import os
import re
import csv
import glob
import shutil
import copy
import json
from datetime import datetime
from taca.utils.filesystem import chdir, control_fastq_filename
from taca.illumina.Runs import Run
from taca.utils import misc
from flowcell_parser.classes import RunParametersParser, SampleSheetParser, RunParser, DemuxSummaryParser


import logging

logger = logging.getLogger(__name__)


class HiSeqX_Run(Run):

    def __init__(self,  run_dir, samplesheet_folders):
        super(HiSeqX_Run, self).__init__( run_dir, samplesheet_folders)
        self._set_sequencer_type()
        self._set_run_type()
        self._copy_samplesheet()

    def _set_sequencer_type(self):
        self.sequencer_type = "HiSeqX"

    def _set_run_type(self):
        self.run_type = "NGI-RUN"

    def _copy_samplesheet(self):
        ssname   = self._get_samplesheet()
        ssparser = SampleSheetParser(ssname)
        try:
            indexfile = self.CONFIG['bcl2fastq']['index_path']
        except KeyError:
            logger.error("Path to index file (10X) not found in the config file")
            raise RuntimeError
        #samplesheet need to be positioned in the FC directory with name SampleSheet.csv (Illumina default)
        #if this is not the case then create it and take special care of modification to be done on the SampleSheet
        samplesheet_dest = os.path.join(self.run_dir, "SampleSheet.csv")
        #Function that goes through the original sample sheet and check for sample types
        self.sample_table = classify_samples(indexfile, ssparser)


                        ################### This part should be removed###############
                        #Function that returns a list of which lanes contains 10X samples.
                        (self.lanes_10X,self.lanes_not_10X,self.tenX_samples) = look_for_lanes_with_10X_indicies(indexfile, ssparser)
                        ##############################################################

        #check that the samplesheet is not already present. In this case go the next step
        if not os.path.exists(samplesheet_dest):
            try:
                with open(samplesheet_dest, 'wb') as fcd:
                    fcd.write(_generate_clean_samplesheet(ssparser,indexfile, rename_samples=True, rename_qPCR_suffix = True, fields_qPCR=[ssparser.dfield_snm]))
            except Exception as e:
                logger.error("encountered the following exception '{}'".format(e))
                return False
            logger.info(("Created SampleSheet.csv for Flowcell {} in {} ".format(self.id, samplesheet_dest)))
        ##SampleSheet.csv generated

        ##when demultiplexing SampleSheet.csv is the one I need to use
        ## Need to rewrite so that SampleSheet_0.csv is always used.
        self.runParserObj.samplesheet  = SampleSheetParser(os.path.join(self.run_dir, "SampleSheet.csv"))
        if not self.runParserObj.obj.get("samplesheet_csv"):
            self.runParserObj.obj["samplesheet_csv"] = self.runParserObj.samplesheet.data

    def demultiplex_run(self):
        """
           Demultiplex a Xten run:
            - Make sub-samplesheet based on sample classes
            - Decide correct bcl2fastq command parameters based on sample classes
            - run bcl2fastq conversion
        """

########################
# Test sample_table
sample_table = {
    '1':[{'101':{'type':'ordinary','length':[8,0]}},{'102':{'type':'ordinary','length':[8,8]}},{'103':{'type':'ordinary','length':[8,8]}},{'104':{'type':'ordinary','length':[8,0]}}],
    '2':[{'105':{'type':'ordinary','length':[8,0]}},{'106':{'type':'ordinary','length':[6,0]}},{'107':{'type':'ordinary','length':[8,0]}},{'108':{'type':'ordinary','length':[8,0]}}],
    '3':[{'105':{'type':'10X','length':[8,0]}},{'106':{'type':'10X','length':[8,0]}},{'107':{'type':'ordinary','length':[8,0]}},{'108':{'type':'ordinary','length':[8,8]}}],
    '4':[{'105':{'type':'10X','length':[8,0]}},{'106':{'type':'ordinary','length':[8,0]}},{'107':{'type':'ordinary','length':[8,0]}},{'108':{'type':'ordinary','length':[8,0]}}]
}
lane_table = {
    '1': [[8,0]],
    '2': [[8,0],[8,8]],
    '3': [[8,0],[8,8],[6,0]],
    '4': [[6,0]],
    '5': [[8,0],[8,8],[6,0],[6,6],[7,0]],
    '6': [[6,0]]
}
########################
        # Check sample types
        sample_type_list = []
        for lane, lane_contents in self.sample_table.items():
            for sample in lane_contents:
                if sample[sample.keys()[0]]['type'] not in sample_type_list:
                    sample_type_list.append(sample[sample.keys()[0]]['type'])

        # Go through sample_table for demultiplexing
        bcl2fastq_cmd_counter = 0
        for type in sample_type_list:
            # Looking for lanes with multiple masks under the same sample type
            lane_table = dict()
            for lane, lane_contents in self.sample_table.items():
                for sample in lane_contents:
                    if sample[sample.keys()[0]]['type'] == type:
                        if lane_table.get(lane):
                            if sample[sample.keys()[0]]['index_length'] not in lane_table[lane]:
                                lane_table[lane].append(sample[sample.keys()[0]]['index_length'])
                        else:
                            lane_table.update({lane:[sample[sample.keys()[0]]['index_length']]})
            # Determine the number of demux needed for the same sample type
            demux_number_with_the_same_sample_type = len(max([v for k, v in lane_table.items()],key=len))
            # Prepare sub-samplesheets, masks and commands
            for i in range(0,demux_number_with_the_same_sample_type):
                bcl2fastq_cmd_counter += i
                # Prepare sub-samplesheet
                # A dictionary with lane and sample IDs to include
                samples_to_include = dict()
                for lane, lane_contents in self.sample_table.items():
                    try:
                        index_length = lane_table[lane][i]
                        for sample in lane_contents:
                            if sample[sample.keys()[0]]['type'] == type and sample[sample.keys()[0]]['index_length'] == index_length:
                                if samples_to_include.get(lane):
                                    samples_to_include[lane].append(sample.keys()[0])
                                else:
                                    samples_to_include.update({lane:[sample.keys()[0]]})
                    except IndexError:
                        continue

                # Make sub-samplesheet
                with chdir(self.run_dir):
                    samplesheet_dest="SampleSheet_{}.csv".format(bcl2fastq_cmd_counter)
                    with open(samplesheet_dest, 'wb') as fcd:
                        fcd.write(_generate_samplesheet_subset(self.runParserObj.samplesheet, samples_to_include))

                # Make base masks
                per_lane_base_masks = self._generate_per_lane_base_mask(type, samples_to_include)

                ###### Need to contiue from here!!  _generate_per_lane_base_mask needs to be modified!!#####












        #we have 10x lane - need to split the  samples sheet and build a 10x command for bcl2fastq
        Complex_run = False
        if len(self.lanes_10X) and len(self.lanes_not_10X):
             Complex_run = True

        if Complex_run:
            with chdir(self.run_dir):
                samplesheet_dest_not_10X="SampleSheet_0.csv"
                with open(samplesheet_dest_not_10X, 'wb') as fcd:
                    fcd.write(_generate_samplesheet_subset(self.runParserObj.samplesheet, self.lanes_not_10X, self.tenX_samples, is_10X = False))
                samplesheet_dest_10X="SampleSheet_1.csv"
                with open(samplesheet_dest_10X, 'wb') as fcd:
                    fcd.write(_generate_samplesheet_subset(self.runParserObj.samplesheet, self.lanes_10X, self.tenX_samples, is_10X = True))
        else:
            with chdir(self.run_dir):
                samplesheet_dest="SampleSheet_0.csv"
                with open(samplesheet_dest, 'wb') as fcd:
                    if len(self.lanes_10X):
                        fcd.write(_generate_samplesheet_subset(self.runParserObj.samplesheet, self.lanes_10X, self.tenX_samples, is_10X = True))
                    else:
                        fcd.write(_generate_samplesheet_subset(self.runParserObj.samplesheet, self.lanes_not_10X, self.tenX_samples, is_10X = False))

        per_lane_base_masks = self._generate_per_lane_base_mask()
        max_different_base_masks =  max([len(per_lane_base_masks[base_masks]) for base_masks in per_lane_base_masks])
        if max_different_base_masks > 1:
            # in a HiSeqX run I cannot have different index sizes in the SAME lane
            logger.error("In FC {} found one or more lane with more than one base mask (i.e., different index sizes in \
                         in the same lane".format(self.id))
            return False
        bcl2fastq_cmd_counter = 0
        with chdir(self.run_dir):
            # create Demultiplexing dir, this changes the status to IN_PROGRESS
            if not os.path.exists("Demultiplexing"):
                os.makedirs("Demultiplexing")
        with chdir(self.run_dir):
            if self.lanes_not_10X:
               cmd_normal = self.generate_bcl_command(self.lanes_not_10X, bcl2fastq_cmd_counter)
               misc.call_external_command_detached(cmd_normal, with_log_files = True, prefix="demux_{}".format(bcl2fastq_cmd_counter))
               logger.info(("BCL to FASTQ conversion and demultiplexing started for "
                   "normal run {} on {}".format(os.path.basename(self.id), datetime.now())))
               bcl2fastq_cmd_counter += 1
            if self.lanes_10X:
               cmd_10X = self.generate_bcl_command(self.lanes_10X, bcl2fastq_cmd_counter, is_10X = True)
               misc.call_external_command_detached(cmd_10X, with_log_files = True, prefix="demux_{}".format(bcl2fastq_cmd_counter))
               logger.info(("BCL to FASTQ conversion and demultiplexing started for "
                   "10X run {} on {}".format(os.path.basename(self.id), datetime.now())))
               bcl2fastq_cmd_counter += 1
        return True


    def _aggregate_demux_results(self):
        """
        Take the Stats.json files from the different demultiplexing folders and merges them into one
        """
        ssname   = self._get_samplesheet()
        ssparser = SampleSheetParser(ssname)
        try:
            indexfile = self.CONFIG['bcl2fastq']['index_path']
        except KeyError:
            logger.error("Path to index file (10X) not found in the config file")
            raise RuntimeError
        # Function that returns a list of which lanes contains 10X samples.
        (lanes_10X,lanes_not_10X,tenX_samples) = look_for_lanes_with_10X_indicies(indexfile, ssparser)
        # Get pure 10X or ordinary lanes plus mixed lanes
        lanes_mixed = list(set(lanes_10X).intersection(set(lanes_not_10X)))

        lanes_10X_dict = {}
        lanes_not_10X_dict = {}

        for lane in lanes_10X:
            lanes_10X_dict[lane] = 0
        for lane in lanes_not_10X:
            lanes_not_10X_dict[lane] = 0

        # When there is no mixed lanes
        if len(lanes_mixed) == 0:
            if len(lanes_not_10X_dict) == 0:
                #in this case I have only 10X lanes, so I can treat it 10X lanes as the easy ones
                self._aggregate_demux_results_simple_complex(lanes_10X_dict, {})
            else:
                self._aggregate_demux_results_simple_complex(lanes_not_10X_dict, lanes_10X_dict)
        # When there is mixed lanes. TACA will deal with lanes
        else:
            self._aggregate_demux_results_simple_complex(lanes_not_10X_dict, lanes_10X_dict)


    def generate_bcl_command(self, lanes, bcl2fastq_cmd_counter, is_10X=False):
        #I have everything to run demultiplexing now.
        logger.info('Building a bcl2fastq command')
        per_lane_base_masks = self._generate_per_lane_base_mask()
        with chdir(self.run_dir):
            cl = [self.CONFIG.get('bcl2fastq')['bin']]
            output_dir = "Demultiplexing_{}".format(bcl2fastq_cmd_counter)
            cl.extend(["--output-dir", output_dir])
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            cl_options = []
            if self.CONFIG.get('bcl2fastq').has_key('options'):
                for option in self.CONFIG['bcl2fastq']['options']:
                    cl_options.extend([option])
                # Add the extra 10X command options if we have a 10X run
                if is_10X:
                    cl_options.extend(self.CONFIG['bcl2fastq']['options_10X'])
                # Append all options that appear in the configuration file to the main command.
                for option in cl_options:
                    if isinstance(option, dict):
                        opt, val = option.items()[0]
                        if "output-dir" not in opt:
                            cl.extend(['--{}'.format(opt), str(val)])
                    else:
                        cl.append('--{}'.format(option))

            cl.extend(["--sample-sheet",  os.path.join(os.path.join(self.run_dir, "SampleSheet_{}.csv".format(bcl2fastq_cmd_counter)))])
            #now add the base_mask for each lane
            for lane in sorted(lanes):
                #Iterate thorugh each lane and add the correct --use-bases-mask for that lane
                base_mask = [per_lane_base_masks[lane][bm]['base_mask'] for bm in per_lane_base_masks[lane]][0] # get the base_mask
                base_mask_expr = "{}:".format(lane) + ",".join(base_mask)
                cl.extend(["--use-bases-mask", base_mask_expr])
        return cl


def _generate_clean_samplesheet(ssparser, indexfile, fields_to_remove=None, rename_samples=True, rename_qPCR_suffix = False, fields_qPCR= None):
    """
        Will generate a 'clean' samplesheet, the given fields will be removed.
        if rename_samples is True, samples prepended with 'Sample_'  are renamed to match the sample name
        Will also replace 10X idicies like SI-GA-A3 with proper indicies like TGTGCGGG
    """
    output=""
    ##expand the ssparser if there are 10X lanes
    index_dict=parse_10X_indexes(indexfile) #read the 10X indices
    # Replace 10X index with the 4 actual indicies.
    for sample in ssparser.data:
        if sample['index'] in index_dict.keys():
            x=0
            while x<3:
                new_sample=dict(sample)
                new_sample['index']=index_dict[sample['index']][x]
                ssparser.data.append(new_sample)
                x+=1
            #Set the original 10X index to the 4th correct index
            sample['index']=index_dict[sample['index']][x]

    #Sort to get the added indicies from 10x in the right place
    ssparser.data.sort()

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
            if rename_samples and ssparser.dfield_sid in field :
                try:
                    if rename_qPCR_suffix and ssparser.dfield_snm in fields_qPCR:
                        #substitute SampleID with SampleName, add Sample_ as prefix and remove __qPCR_ suffix
                        value =re.sub('__qPCR_$', '', 'Sample_{}'.format(line[ssparser.dfield_snm]))
                    else:
                        #substitute SampleID with SampleName, add Sample_ as prefix
                        value ='Sample_{}'.format(line[ssparser.dfield_snm])
                except:
                        #otherwise add Sample_ as prefix
                        value = 'Sample_{}'.format(line[ssparser.dfield_sid])
            elif rename_qPCR_suffix and field in fields_qPCR:
                value = re.sub('__qPCR_$', '', line[field])

            line_ar.append(value)

        output+=",".join(line_ar)
        output+=os.linesep

    return output


def classify_samples(indexfile, ssparser):
    """
    Given an ssparser object
    goes through all samples and decide sample types
    """
    TENX_GENO_PAT = re.compile("SI-GA-[A-H][1-9][0-2]?")
    TENX_ATAC_PAT = re.compile("SI-NA-[A-H][1-9][0-2]?")
    UMI_IDX = re.compile("([ATCG]{4,}N+$)")
    sample_table = dict()
    index_dict = parse_10X_indexes(indexfile)
    for sample in ssparser.data:
        lane = sample['Lane']
        sample_name = sample.get('Sample_Name') or sample.get('SampleName')
        # 10X Genomic DNA & RNA
        if TENX_GENO_PAT.findall(sample['index']):
            index_length = [len(index_dict[sample['index']][0]),0]
            sample_type = '10X_GENO'
        # 10X scATAC
        elif TENX_ATAC_PAT.findall(sample['index']):
            index_length = [len(index_dict[sample['index']][0]),16]
            sample_type = '10X_ATAC'
        # UMI samples
        elif UMI_IDX.findall(sample['index']) or UMI_IDX.findall(sample['index2']):
            # Index length after removing "N" part
            index_length = [len(sample['index'].replace('N','')),len(sample['index2'].replace('N',''))]
            sample_type = 'UMI'
        # No Index case. Note that if both index 1 and 2 are empty, it will be the same index type but will be handled in the next case
        elif sample['index'].upper() == 'NOINDEX':
            index_length = [0,0]
            sample_type = 'ordinary'
        # Ordinary samples
        else:
            index_length = [len(sample['index']),len(sample['index2'])]
            sample_type = 'ordinary'

        # Write in sample table
        if sample_table.get(lane):
            sample_table[lane].append({sample_name:{'type':sample_type,'index_length':index_length}})
        else:
            sample_table.update({lane:[{sample_name:{'type':sample_type,'index_length':index_length}}]})

        return sample_table


def look_for_lanes_with_10X_indicies(indexfile, ssparser):
    """
    Given an ssparser object
    returns a list of lanes with 10X indicies
    """
    index_dict=parse_10X_indexes(indexfile)
    tenX_lanes = set() #Set to only get each lane once
    not_tenX_lanes = set()
    tenX_samples = dict()
    for sample in ssparser.data:
        lane = sample['Lane']
        if sample['index'] in index_dict.keys():
            tenX_lanes.add(lane)
            sample_name = sample.get('Sample_Name') or sample.get('SampleName')
            if tenX_samples.get(lane):
                tenX_samples[lane].append(sample_name)
            else:
                tenX_samples.update({lane:[sample_name]})
        else:
            not_tenX_lanes.add(lane)

    return (list(tenX_lanes), list(not_tenX_lanes), tenX_samples)


def parse_10X_indexes(indexfile):
    """
    Takes a file of 10X indexes and returns them as a dict.
    Todo: Set it up to take the file from config instead
    """
    index_dict={}
    with open(indexfile , 'r') as f:
        for line in f:
            line_=line.rstrip().split(',')
            index_dict[line_[0]]=line_[1:5]
    return index_dict


def _generate_samplesheet_subset(ssparser, samples_to_include):
    output=""

    #Header
    output+="[Header]{}".format(os.linesep)
    for field in ssparser.header:
        output+="{},{}".format(field.rstrip(), ssparser.header[field].rstrip())
        output+=os.linesep
    #Data
    output+="[Data]{}".format(os.linesep)
    datafields=[]
    for field in ssparser.datafields:
        datafields.append(field)
    output+=",".join(datafields)
    output+=os.linesep
    for line in ssparser.data:
        sample_name = line.get('Sample_Name') or line.get('SampleName')
        lane = line['Lane']
        if lane in samples_to_include.keys():
            if sample_name in samples_to_include.get(lane):
                line_ar=[]
                for field in datafields:
                    # Case of no index
                    if field == "index" and "NOINDEX" in line[field].upper():
                        line[field] = ""
                    # Case of UMI
                    if (field == "index" or field == "index2") and UMI_IDX.findall(line[field]):
                        line[field] = line[field].replace('N','')
                    line_ar.append(line[field])
                output+=",".join(line_ar)
                output+=os.linesep

    return output


def _generate_per_lane_base_mask(self):
    """
    This functions generate the base mask for each lane.
    Hypotesis:
        - RunInfo.xml contains the configuration
        - this object contains a properly parsed samplesheet
    It returns an dict with a key for each lane:
    {lane1:
        {base_mask_string (e.g., Y150I6N2N8Y150):
            [ base_mask , [SampleSheetEntries]]
        }
     lane2:
    }
    """
    # generate new ssparser (from the renamed smaplesheet)
    runSetup = self.runParserObj.runinfo.get_read_configuration()
    base_masks = {}
    if not self.runParserObj.samplesheet:
        raise RuntimeError("samplesheet not yet initialised")

    for data_entry in self.runParserObj.samplesheet.data:
        ## for each data_entry in my samplesheet (i.e., for each sample)
        lane  = data_entry['Lane']
        if lane not in base_masks:
            base_masks[lane] = {}
        index = ""
        index2 = ""
        is_dual_index = False
        if data_entry.get('index'):
            index = data_entry['index']
            if index in "NoIndex": #special case for HiSeq when one sample is alone in a lane
                index = ""
            is_dual_index = False # default for Xten
            if data_entry.get('index2'):
                index2 = data_entry['index2']
                is_dual_index = True
            #specific for HiSeq, will disapper once we will use bcl2fastq_2.17
            #index = data_entry['Index'].replace('-', '').replace('NoIndex', '')
        index_size  = len(index)
        index2_size = len(index2)
        # compute the basemask
        base_mask = self._compute_base_mask(runSetup, index_size, is_dual_index, index2_size)
        base_mask_string = "".join(base_mask)
        # prepare the dictionary
        if base_mask_string not in base_masks[lane]:
            # first time I find such base mask in this lane,
            base_masks[lane][base_mask_string] = {'base_mask':base_mask,
                                                  'data' : []}
        base_masks[lane][base_mask_string]['data'].append(data_entry)

    return base_masks


def _compute_base_mask(self, runSetup, index_size, dual_index_sample, index2_size):
    """
        Assumptions:
            - if runSetup is of size 3, then single index run
            - if runSetup is of size 4, then dual index run
    """
    bm = []
    dual_index_run = False
    if len(runSetup) > 4:
        raise RuntimeError("when generating base_masks looks like there are"
                           " more than 4 reads in the RunSetup.xml")

    for read in runSetup:
        cycles = int(read['NumCycles'])
        if read['IsIndexedRead'] == 'N':
            bm.append('Y' + str(cycles))
        else:
            if index_size > cycles:
                # the size of the index of the sample sheet is larger than the
                # one specified by RunInfo.xml, somethig must be wrong
                raise RuntimeError("when generating base_masks found index in"
                                   "samplesheet larger than the index specifed in RunInfo.xml")
            is_first_index_read = int(read['Number']) == 2
            # now prepare the base mask for this index read
            if is_first_index_read:
                i_remainder = cycles - index_size
                if i_remainder > 0:
                    if index_size == 0:
                        bm.append('N' + str(cycles)) #special case (NoIndex)
                    else:
                        bm.append('I' + str(index_size) + 'N' + str(i_remainder))
                else:
                    bm.append('I' + str(cycles))
            else:
            # when working on the second read index I need to know if the sample is dual index or not
                if dual_index_sample:
                    i_remainder = cycles - index2_size
                    if i_remainder > 0:
                        if index2_size == 0:
                            bm.append('N' + str(cycles)) #possible if same lane has single and dual index samples
                        else:
                            bm.append('I' + str(index2_size) + 'N' + str(i_remainder))
                    else:
                        bm.append('I' + str(cycles))
                else:
                # if this sample is not dual index but the run is,
                # then I need to ignore the second index completely
                    bm.append('N' + str(cycles))
    return bm
