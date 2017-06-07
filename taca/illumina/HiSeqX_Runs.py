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



    def _set_sequencer_type(self):
        self.sequencer_type = "HiSeqX"

    def _set_run_type(self):
        self.run_type = "NGI-RUN"


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
        ##DEBUING - ReMOVE
        self.aggregate_results()
        return 
        ##END of DEBUGGING
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
        #Function that returns a list of which lanes contains 10X samples.
        (lanes_10X,lanes_not_10X) = look_for_lanes_with_10X_indicies(indexfile, ssparser)
        #check that the samplesheet is not already present. In this case go the next step
        if not os.path.exists(samplesheet_dest):
            try:
                with open(samplesheet_dest, 'wb') as fcd:
                    fcd.write(_generate_clean_samplesheet(ssparser,indexfile, fields_to_remove=['index2'], rename_samples=True, rename_qPCR_suffix = True, fields_qPCR=[ssparser.dfield_snm]))
            except Exception as e:
                logger.error("encountered the following exception '{}'".format(e))
                return False
            logger.info(("Created SampleSheet.csv for Flowcell {} in {} ".format(self.id, samplesheet_dest)))
        ##SampleSheet.csv generated

        ##when demultiplexing SampleSheet.csv is the one I need to use
        ## Need to rewrite so that SampleSheet_0.csv is always used.
        self.runParserObj.samplesheet  = SampleSheetParser(os.path.join(self.run_dir, "SampleSheet.csv"))
        #we have 10x lane - need to split the  samples sheet and build a 10x command for bcl2fastq
        Complex_run = False
        if len(lanes_10X) and len(lanes_not_10X):
            Complex_run = True

        if Complex_run:
            samplesheet_dest_not_10X="SampleSheet_0.csv"
            with open(samplesheet_dest_not_10X, 'wb') as fcd:
                fcd.write(_generate_samplesheet_subset(self.runParserObj.samplesheet, lanes_not_10X))
            samplesheet_dest_10X="SampleSheet_1.csv"
            with open(samplesheet_dest_10X, 'wb') as fcd:
                fcd.write(_generate_samplesheet_subset(self.runParserObj.samplesheet, lanes_10X))
        else:
            with chdir(self.run_dir):
                shutil.copy("SampleSheet.csv", "SampleSheet_0.csv")
        
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
            if lanes_not_10X:
               cmd_normal = self.generate_bcl_command(lanes_not_10X, bcl2fastq_cmd_counter)
               misc.call_external_command_detached(cmd_normal, with_log_files = True, prefix="demux_{}".format(bcl2fastq_cmd_counter))
               logger.info(("BCL to FASTQ conversion and demultiplexing started for "
                   "normal run {} on {}".format(os.path.basename(self.id), datetime.now())))
               bcl2fastq_cmd_counter += 1
            if lanes_10X:
               cmd_10X = self.generate_bcl_command(lanes_10X, bcl2fastq_cmd_counter, is_10X = True)
               misc.call_external_command_detached(cmd_10X, with_log_files = True, prefix="demux_{}".format(bcl2fastq_cmd_counter))
               logger.info(("BCL to FASTQ conversion and demultiplexing started for "
                   "10X run {} on {}".format(os.path.basename(self.id), datetime.now())))
               bcl2fastq_cmd_counter += 1
        return True

    
    def aggregate_results(self):
        """
        Take the Stats.json files from the different demultiplexing folders and merges them into one
        """
        complex_run = False 
        os.chdir(self.run_dir)
        if os.path.exists("Demultiplexing_1") and os.path.exists("Demultiplexing_0"):
            complex_run = True 
            # Results needs to be aggregated 
        elif os.path.exists("Demultiplexing_0"):
            #Simple run, no need to aggregate results.
            #But Demultiplexing_0 needs to be symlinked into just Demultiplexing
            os.rmdir("Demultiplexing")
            os.symlink("Demultiplexing_0", "Demultiplexing") 
        else:
            logger.error("Could not open Demultiplexing directory")
            return
        if complex_run:
            if not os.path.exists("Demultiplexing"):
                os.makedirs("Demultiplexing")  ##Making it if it was removed, only useful while debugging

            ## We need to merge the Stats.json files and the Demltiplexing directories 
            stats_list=[] #List with the dicts of the two 
            dir_num=0
            while dir_num < 2:
                os.chdir("Demultiplexing_{}/Stats".format(dir_num))
                with open('Stats.json') as json_data:
                     data = json.load(json_data)
                     stats_list.append(data)
                dir_num +=1
                os.chdir(self.run_dir)
            #Update the info in Demux_0 (Normal) with the info in Demux_1 (10X)
            stats_list[0]['ReadInfosForLanes'].extend(stats_list[1]['ReadInfosForLanes'])
            stats_list[0]['ConversionResults'].extend(stats_list[1]['ConversionResults'])
            stats_list[0]['UnknownBarcodes'].extend(stats_list[1]['UnknownBarcodes'])
            stats_list.sort()
            #Now we need to merge the Demultiplexing_ dirs into Demultiplexing. 
            demux_dirs = glob.glob(os.path.join(self.run_dir,"Demultiplexing_*")) 
            dest = os.path.join(self.run_dir, "Demultiplexing")
            src_files=[]
            for sub_dir in demux_dirs:
                 src_files.extend(glob.glob(os.path.join(demux_dirs[0],"*")))
            for file_name in src_files:
                if 'Reports' in file_name  or 'Stats' in file_name:        ## The reports are specific to each run and would overwrite each other.            
                    continue 
                if (os.path.isdir(file_name)) or os.path.isfile(file_name):
                    if  os.path.basename(os.path.dirname(file_name)).startswith("Demultiplexing"):
                        try:
                            os.symlink(file_name, os.path.join(dest, os.path.basename(file_name)))
                        except Exception as e:
                             if e.errno == 17:
                                 continue  #Don't want an error is file exists, otherwise this will fire for all undetrmined files every time
                             else:
                                 logger.info("While trying to create a symlink for {}, TACA encountered error: '{}'".format(file_name,e))
                    elif os.path.isfile(file_name):   
                         ## add the subdir of the file to the new path name. i.e Stats to Stats.json. 
                         try:
                            os.symlink(file_name, os.path.join(os.path.join(dest, os.path.basename(os.path.normpath(os.path.dirname(file_name))), os.path.basename(file_name))))
                         except Exception as e:
                             if e.errno == 17:
                                 continue  #Don't want an error is file exists, otherwise this will fire for all undetrmined files every time
                             else:
                                     logger.info("While trying to create a symlink for {}, TACA encountered error: '{}'".format(file_name,e))
            Statsdir = (os.path.join(dest, "Stats"))
            if not os.path.isdir(Statsdir):
                os.mkdir(Statsdir)
            with open (os.path.join(Statsdir, "Stats.json"), 'w') as json_out:
                 json.dump(stats_list[0], json_out)
        
        else:
            return 
             
        
        
    def check_run_status(self):
        """
           This function checks the status of a run while in progress.
            When the run is completed kick of the results aggregation step.
        """
        run_dir    =  self.run_dir
        dex_status =  self.get_run_status()
        return None
        
        ##Copied from Hiseq class
        if  dex_status == 'COMPLETED':
            return None
        #otherwise check the status of running demux
        #collect all samplesheets generated before
        samplesheets =  glob.glob(os.path.join(run_dir, "*_[0-9].csv")) # a single digit... this hipotesis should hold for a while
        allDemuxDone = True
        for samplesheet in samplesheets:
            #fetch the id of this demux job
            demux_id = os.path.splitext(os.path.split(samplesheet)[1])[0].split("_")[1]
            #demux folder is
            demux_folder = os.path.join(run_dir, "Demultiplexing_{}".format(demux_id))
            #check if this job is done
            if os.path.exists(os.path.join(run_dir, demux_folder, 'Stats', 'DemultiplexingStats.xml')):
                allDemuxDone = allDemuxDone and True
                logger.info("Sub-Demultiplexing in {} completed.".format(demux_folder))
            else:
                allDemuxDone = allDemuxDone and False
                logger.info("Sub-Demultiplexing in {} not completed yet.".format(demux_folder))
        #in this case, I need to aggreate in the Demultiplexing folder all the results
        if allDemuxDone:
            self.aggregate_results()
            #now I can initialise the RunParser
            self.runParserObj = RunParser(self.run_dir)


    def compute_undetermined(self):
        """
            This function parses the Undetermined files per lane produced by illumina
            for now nothign done, TODO: check all undetermined files are present as sanity check
        """
        return True

    
    def check_undetermined_reads(self, lane, freq_tresh):
        """checks that the number of undetermined reads does not exceed a given threshold
        returns true if the percentage is lower then freq_tresh
        Does this by considering undetermined all reads marked as unknown

        :param lane: lane identifier
        :type lane: string
        :param freq_tresh: maximal allowed percentage of undetermined indexes in a lane
        :type frew_tresh: float
        :rtype: boolean
        :returns: True if the checks passes, False otherwise
        """
        #compute lane yield
        run = self.run_dir
        lanes = self.runParserObj.lanes
        lane_yield = 0;
        for entry in lanes.sample_data:
            if lane == entry['Lane']:
                if lane_yield > 0:
                    logger.warn("lane_yeld must be 0, somehting wrong is going on here")
                lane_yield = int(entry['PF Clusters'].replace(',',''))

        #I do not need to parse undetermined here, I can use the the lanes object to fetch unknown
        sample_lanes = self.runParserObj.lanebarcodes
        undetermined_lane_stats = [item for item in sample_lanes.sample_data if item["Lane"]==lane and item["Sample"]=="Undetermined"]
        undetermined_total = 0
        percentage_und     = 0
        if len(undetermined_lane_stats) > 1:
            logger.error("Something wrong in check_undetermined_reads, found more than one undetermined sample in one lane")
            return False
        elif len(undetermined_lane_stats) == 0:
            #NoIndex case
            undetermined_total = 0
            percentage_und     = 0
        else:
            #normal case
            undetermined_total = int(undetermined_lane_stats[0]['PF Clusters'].replace(',',''))
            percentage_und = (undetermined_total/float(lane_yield))*100

        if  percentage_und > freq_tresh:
            logger.warn("The undetermined indexes account for {}% of lane {}, "
                        "which is over the threshold of {}%".format(percentage_und, lane, freq_tresh))
            return False
        else:
            return True


    def check_maximum_undertemined_freq(self, lane, freq_tresh):
        """returns true if the most represented index accounts for less than freq_tresh
            of the total amount of undetermiend

            :param lane: lane identifier
            :type lane: string
            :param freq_tresh: maximal allowed frequency of the most frequent undetermined index
            :type frew_tresh: float
            :rtype: boolean
            :returns: True if the checks passes, False otherwise
            """

        #check the most reptresented index
        undeterminedStats = DemuxSummaryParser(os.path.join(self.run_dir,self.demux_dir, "Stats"))
        most_frequent_undet_index_count = int(undeterminedStats.result[lane].items()[0][1])
        most_frequent_undet_index       = undeterminedStats.result[lane].items()[0][0]

        #compute the total amount of undetermined reads
        sample_lanes = self.runParserObj.lanebarcodes
        undetermined_lane_stats = [item for item in sample_lanes.sample_data if item["Lane"]==lane and item["Sample"]=="Undetermined"]
        freq_most_occuring_undet_index = 0
        if len(undetermined_lane_stats) > 1:
            logger.error("Something wrong in check_undetermined_reads, found more than one undetermined sample in one lane")
            return False
        elif len(undetermined_lane_stats) == 0:
            #NoIndex case
            freq_most_occuring_undet_index = 0
        else:
            undetermined_total = int(undetermined_lane_stats[0]['PF Clusters'].replace(',',''))
            freq_most_occuring_undet_index = (most_frequent_undet_index_count/float(undetermined_total))*100

        if freq_most_occuring_undet_index > freq_tresh:
            logger.warn("The most frequent barcode of lane {} ({}) represents {}%, "
                        "which is over the threshold of {}%".format(lane, most_frequent_undet_index, freq_most_occuring_undet_index , freq_tresh))
            return False
        else:
            return True





    def get_path_per_lane(self):
        """
        :param run: the path to the flowcell
        :type run: str
        :param ss: SampleSheet reader
        :type ss: flowcell_parser.XTenSampleSheet
        """
        run          = self.run_dir
        dmux_folder  = self.demux_dir
        ss           = self.runParserObj.samplesheet
        d={}
        for l in ss.data:
            try:
                d[l['Lane']]=os.path.join(run, dmux_folder, l[ss.dfield_proj], l[ss.dfield_sid])
            except KeyError:
                logger.error("Can't find the path to the sample, is 'Project' in the samplesheet ?")
                d[l['Lane']]=os.path.join(run, dmux_folder)

        return d

    def get_samples_per_lane(self):
        """
        :param ss: SampleSheet reader
        :type ss: flowcell_parser.XTenSampleSheet
        :rtype: dict
        :returns: dictionnary of lane:samplename
        """
        ss = self.runParserObj.samplesheet
        d={}
        for l in ss.data:
            s=l[ss.dfield_snm].replace("Sample_", "").replace("-", "_")
            d[l['Lane']]=l[ss.dfield_snm]

        return d



    def _rename_undet(self, lane, samples_per_lane):
        """Renames the Undetermined fastq file by prepending the sample name in front of it

        :param run: the path to the run folder
        :type run: str
        :param status: the demultiplexing status
        :type status: str
        :param samples_per_lane: lane:sample dict
        :type status: dict
        """
        run = self.run_dir
        dmux_folder = self.demux_dir
        for file in glob.glob(os.path.join(run, dmux_folder, "Undetermined*L0?{}*".format(lane))):
            old_name=os.path.basename(file)
            old_name_comps=old_name.split("_")
            old_name_comps[1]=old_name_comps[0]# replace S0 with Undetermined
            old_name_comps[0]=samples_per_lane[lane]#replace Undetermined with samplename
            for index, comp in enumerate(old_name_comps):
                if comp.startswith('L00'):
                    old_name_comps[index]=comp.replace('L00','L01')#adds a 1 as the second lane number in order to differentiate undetermined from normal in piper

            new_name="_".join(old_name_comps)
            logger.info("Renaming {} to {}".format(file, os.path.join(os.path.dirname(file), new_name)))
            os.rename(file, os.path.join(os.path.dirname(file), new_name))

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
            if self.CONFIG.get('bcl2fastq').has_key('options'):
                cl_options = self.CONFIG['bcl2fastq']['options']
            
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

def look_for_lanes_with_10X_indicies(indexfile, ssparser):
    """
    Given an ssparser object
    returns a list of lanes with 10X indicies
    """
    index_dict=parse_10X_indexes(indexfile)
    tenX_lanes = set() #Set to only get each lane once
    not_tenX_lanes = set()
    for sample in ssparser.data:
        if sample['index'] in index_dict.keys():
            tenX_lanes.add(sample['Lane'])
        else:
            not_tenX_lanes.add(sample['Lane'])

    return (list(tenX_lanes),list(not_tenX_lanes))


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


def _generate_samplesheet_subset(ssparser, lanes):
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
        if line['Lane'] in lanes:
            line_ar=[]
            for field in datafields:
                line_ar.append(line[field])
            output+=",".join(line_ar)
            output+=os.linesep

    return output
