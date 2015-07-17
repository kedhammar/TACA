import os
import re
import csv
import glob
from datetime import datetime
from taca.utils.filesystem import chdir, control_fastq_filename
from taca.illumina.Runs import Run
from taca.utils import misc
from flowcell_parser.classes import RunParametersParser, SampleSheetParser, RunParser
import taca.illumina.HiSeqX_QC as qc


import logging

logger = logging.getLogger(__name__)

class HiSeqX_Run(Run):
    
    def __init__(self,  run_dir, samplesheet_folders):
        super(HiSeqX_Run, self).__init__( run_dir, samplesheet_folders)
        self._set_sequencer_type()


    def _set_sequencer_type(self):
        self.sequencer_type = "HiSeqX"



    def _get_samplesheet(self):
        """
            Locate and parse the samplesheet for a run. The idea is that there is a folder in
            samplesheet_folders that contains a samplesheet named flowecell_id.csv.
        """
        current_year = '20' + self.id[0:2]
        samplesheets_dir = os.path.join(self.CONFIG['samplesheets_dir'],
                                                current_year)
        ssname = os.path.join(samplesheets_dir, '{}.csv'.format(self.flowcell_id))
        if os.path.exists(ssname):
            return ssname
        else:
            raise RuntimeError("not able to find samplesheet {}.csv in {}".format(self.flowcell_id, self.CONFIG['samplesheets_dir']))





    def demultiplex_run(self):
        """
           Demultiplex a Xten run:
            - find the samplesheet
            - make a local copy of the samplesheet and name it SampleSheet.csv
            - define if necessary the bcl2fastq commands (if indexes are not of size 8, i.e. neoprep)
            - run bcl2fastq conversion
        """

        ssname   = self._get_samplesheet()
        ssparser = SampleSheetParser(ssname)


        #samplesheet need to be positioned in the FC directory with name SampleSheet.csv (Illumina default)
        #if this is not the case then create it and take special care of modification to be done on the SampleSheet
        samplesheet_dest = os.path.join(self.run_dir, "SampleSheet.csv")
        #check that the samplesheet is not already present and generate the dafault SampleSheet. This avoids multiple runs on the same FC
        if not os.path.exists(samplesheet_dest):
            try:
                with open(samplesheet_dest, 'wb') as fcd:
                    fcd.write(_generate_clean_samplesheet(ssparser, fields_to_remove=['index2'], rename_samples=True, rename_qPCR_suffix = True, fields_qPCR=['SampleName']))
            except Exception as e:
                logger.error(e.text)
                return False
            logger.info(("Created SampleSheet.csv for Flowcell {} in {} ".format(self.id, samplesheet_dest)))
        ##SampleSheet.csv generated
        ##when demultiplexing SampleSheet.csv is the one I need to use
        self.runParserObj.samplesheet  = SampleSheetParser(os.path.join(self.run_dir, "SampleSheet.csv"))

        per_lane_base_masks = self._generate_per_lane_base_mask()
        max_different_base_masks =  max([len(base_masks) for base_masks in per_lane_base_masks])
        if max_different_base_masks > 1:
            # in a HiSeqX run I cannot have different index sizes in the SAME lane
            logger.error("In FC {} found one or more lane with more than one base mask (i.e., different index sizes in \
                         in the same lane".format(self.id))
            return False
        #I have everything to run demultiplexing now.
        logger.info('Building bcl2fastq command')

        with chdir(self.run_dir):
            cl = [self.CONFIG.get('bcl2fastq')['bin']]
            if self.CONFIG.get('bcl2fastq').has_key('options'):
                cl_options = self.CONFIG['bcl2fastq']['options']
                # Append all options that appear in the configuration file to the main command.
                for option in cl_options:
                    if isinstance(option, dict):
                        opt, val = option.popitem()
                        cl.extend(['--{}'.format(opt), str(val)])
                    else:
                        cl.append('--{}'.format(option))
            #now add the base_mask for each lane
            for lane in sorted(per_lane_base_masks):
                #iterate thorugh each lane and add the correct --use-bases-mask for that lane
                #there is a single basemask for each lane, I checked it a couple of lines above
                base_mask = [per_lane_base_masks[lane][bm]['base_mask'] for bm in per_lane_base_masks[lane]][0] # get the base_mask
                base_mask_expr = "{}:".format(lane) + ",".join(base_mask)
                cl.extend(["--use-bases-mask", base_mask_expr])

            logger.info(("BCL to FASTQ conversion and demultiplexing started for "
                 " run {} on {}".format(os.path.basename(self.id), datetime.now())))
            misc.call_external_command_detached(cl, with_log_files=True)
        return True


    def check_run_status(self):
        """
           This function checks the status of a run while in progress.
           For Xten analysis it also compute undetermined stats
        """
        run_dir    =  self.run_dir
        dex_status =  self.get_run_status()
        qc.compute_undetermined_stats(self.run_dir, self.demux_dir, self.get_run_status())

    def demux_done(self):
        """
           checks that the demux is done and that there are no concurrent TACA processes running
        """
        if self.get_run_status() is not 'COMPLETED':
            return False
        #if demux is completed check that there are no concurrent executions of TACA
        running = 0
        for file in glob.glob(os.path.join(self.run_dir, self.demux_dir, "index_count_L*.running")):
            running +=1
        if running > 0:
            return False
        #replace - with _ between project and sample name in the fastq
        control_fastq_filename(os.path.join(self.run_dir, self.demux_dir))
        return True

    def post_demux(self):
        raise NotImplementedError("Please Implement this method")



    def check_QC(self):
        if self.demux_done():
            return qc.check_lanes_QC(self)
        else:
            raise RuntimeError("Trying to QC a run but the run is not compelted yet")








def _generate_clean_samplesheet(ssparser, fields_to_remove=None, rename_samples=True, rename_qPCR_suffix = False, fields_qPCR= None):
    """
        Will generate a 'clean' samplesheet, the given fields will be removed. 
        if rename_samples is True, samples prepended with 'Sample_'  are renamed to match the sample name
    """
    output=""
    if not fields_to_remove:
        fields_to_remove=[]
    #Header
    output+="[Header]{}".format(os.linesep)
    for field in ssparser.header:
        output+="{},{}".format(field.rstrip(), ssparser.header[field].rstrip())
        output+=os.linesep
    #Data
    output+="[Data]{}".format(os.linesep)
    datafields=[]
    for field in ssparser.datafields:
        if field not in fields_to_remove:
            datafields.append(field)
    output+=",".join(datafields)
    output+=os.linesep
    for line in ssparser.data:
        line_ar=[]
        for field in datafields:
            value = line[field]
            if rename_samples and 'SampleID' in field :
                try:
                    if rename_qPCR_suffix and 'SampleName' in fields_qPCR:
                        #substitute SampleID with SampleName, add Sample_ as prefix and remove __qPCR_ suffix
                        value =re.sub('__qPCR_$', '', 'Sample_{}'.format(line['SampleName']))
                    else:
                        #substitute SampleID with SampleName, add Sample_ as prefix
                        value ='Sample_{}'.format(line['SampleName'])
                except:
                        #otherwise add Sample_ as prefix
                        value = 'Sample_{}'.format(line['SampleID'])
            elif rename_qPCR_suffix and field in fields_qPCR:
                value = re.sub('__qPCR_$', '', line[field])
                                                                                                                            
            line_ar.append(value)
                                                                                                                                
        output+=",".join(line_ar)
        output+=os.linesep

    return output

    