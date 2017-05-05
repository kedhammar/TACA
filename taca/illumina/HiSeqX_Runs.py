import os
import re
import csv
import glob
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

        ssname   = self._get_samplesheet()
        ssparser = SampleSheetParser(ssname)
        #samplesheet need to be positioned in the FC directory with name SampleSheet.csv (Illumina default)
        #if this is not the case then create it and take special care of modification to be done on the SampleSheet
        samplesheet_dest = os.path.join(self.run_dir, "SampleSheet.csv")
        #Function that returns a list of which lanes contains 10X samples.
        (lanes_10X,lanes_not_10X) = look_for_lanes_with_10X_indicies(ssparser)
        #check that the samplesheet is not already present. In this case go the next step
        if not os.path.exists(samplesheet_dest):
            try:
                with open(samplesheet_dest, 'wb') as fcd:
                    fcd.write(_generate_clean_samplesheet(ssparser, fields_to_remove=['index2'], rename_samples=True, rename_qPCR_suffix = True, fields_qPCR=[ssparser.dfield_snm]))
            except Exception as e:
                logger.error(e.text)
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
            shutil.copyfile("SampleSheet.csv", "SampleSheet_0.csv")

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
        if lanes_not_10X:
            cmd_normal = self.generate_bcl_command(lanes_not_10X, bcl2fastq_cmd_counter)
            misc.call_external_command_detached(cmd_normal, with_log_files = True )
            logger.info(("BCL to FASTQ conversion and demultiplexing started for "
                "normal run {} on {}".format(os.path.basename(self.id), datetime.now())))
            bcl2fastq_cmd_counter += 1

        if lanes_10X:
            cmd_10X = self.generate_bcl_command(lanes_10X, bcl2fastq_cmd_counter, is_10X = True)
            misc.call_external_command_detached(cmd_10X, with_log_files = True )
            logger.info(("BCL to FASTQ conversion and demultiplexing started for "
                " 10X run {} on {}".format(os.path.basename(self.id), datetime.now())))
            bcl2fastq_cmd_counter += 1

        return True


    def check_run_status(self):
        """
           This function checks the status of a run while in progress. In the future will print the status
        """
        run_dir    =  self.run_dir
        dex_status =  self.get_run_status()
        return None


    def compute_undetermined(self):
        """
            This function parses the Undetermined files per lane produced by illumina
            for now nothign done, TODO: check all undetermined files are present as sanity check
        """
        return True



    def check_QC(self):
        run_dir = self.run_dir
        dmux_folder = self.demux_dir

        max_percentage_undetermined_indexes_pooled_lane   = self.CONFIG['QC']['max_percentage_undetermined_indexes_pooled_lane']
        max_percentage_undetermined_indexes_unpooled_lane = self.CONFIG['QC']['max_percentage_undetermined_indexes_unpooled_lane']
        minimum_percentage_Q30_bases_per_lane             = self.CONFIG['QC']['minimum_percentage_Q30_bases_per_lane']
        minimum_yield_per_lane                            = self.CONFIG['QC']['minimum_yield_per_lane']
        max_frequency_most_represented_und_index_pooled_lane   = self.CONFIG['QC']['max_frequency_most_represented_und_index_pooled_lane']
        max_frequency_most_represented_und_index_unpooled_lane = self.CONFIG['QC']['max_frequency_most_represented_und_index_unpooled_lane']

        if not self.runParserObj.samplesheet or not self.runParserObj.lanebarcodes or not self.runParserObj.lanes:
            logger.error("Something went wrong while parsing demultiplex results. QC cannot be performed.")
            return False

        status = True #initialise status as passed
        #read the samplesheet and fetch all lanes
        lanes_to_qc       = misc.return_unique([lanes['Lane'] for lanes in  self.runParserObj.samplesheet.data])
        path_per_lane    =  self.get_path_per_lane()
        samples_per_lane =  self.get_samples_per_lane()
        #now for each lane
        for lane in lanes_to_qc:
            lane_status = True
            #QC lane yield
            if self.lane_check_yield(lane, minimum_yield_per_lane):
                lane_status = lane_status and True
            else:
                logger.warn("lane {} did not pass yield qc check. This FC will not be transferred.".format(lane))
                lane_status = lane_status and False
            #QC on the total %>Q30 of the all lane
            if self.lane_check_Q30(lane, minimum_percentage_Q30_bases_per_lane):
                lane_status = lane_status and True
            else:
                logger.warn("lane {} did not pass Q30 qc check. This FC will not be transferred.".format(lane))
                lane_status = lane_status and False
            #QC for undetermined
            max_percentage_undetermined_indexes = max_percentage_undetermined_indexes_pooled_lane
            max_frequency_most_represented_und  = max_frequency_most_represented_und_index_pooled_lane
            #distinguish the case between Pooled and Unpooled lanes, for unpooled lanes rename the Undetemriend file
            if self.is_unpooled_lane(lane):
                ##DO NOT ADD UNDET BY DEFAULT TO SAMPLES
                #rename undetermiend, in this way PIPER will be able to use them
                self._rename_undet(lane, samples_per_lane)
                ##logger.info("linking undetermined lane {} to sample".format(lane))
                ##but do not soft link them
                #misc.link_undet_to_sample(run_dir, dmux_folder, lane, path_per_lane)
                max_percentage_undetermined_indexes = max_percentage_undetermined_indexes_unpooled_lane
                max_frequency_most_represented_und  = max_frequency_most_represented_und_index_unpooled_lane


            if self.check_undetermined_reads(lane, max_percentage_undetermined_indexes):
                if self.check_maximum_undertemined_freq(lane, max_frequency_most_represented_und):
                    lane_status= lane_status and True
                else:
                    logger.warn("lane {} did not pass the check for most represented undet index. Most occuring undetermined index occurs too often.".format(lane))
                    lane_status= lane_status and False
            else:
                logger.warn("lane {} did not pass the undetermined qc checks. Fraction of undetermined too large.".format(lane))
                lane_status= lane_status and False
            if lane_status:
                logger.info("lane {} passed all qc checks".format(lane))
            #store the status for the all FC
            status = status and lane_status

        return status




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
        logger.info('Building bcl2fastq command for normal lanes')
        per_lane_base_masks = self._generate_per_lane_base_mask()
        with chdir(self.run_dir):
            if not os.path.exists("Demultiplexing_{}".format(bcl2fastq_cmd_counter)):
                os.makedirs("Demultiplexing_{}".format(bcl2fastq_cmd_counter))
            cl = [self.CONFIG.get('bcl2fastq')['bin']]
            if self.CONFIG.get('bcl2fastq').has_key('options'):
                cl_options = self.CONFIG['bcl2fastq']['options']
                # Add the extra 10X command options if we have a 10X run
                if is_10X:
                    cl_options.append(self.CONFIG['bcl2fastq']['options_10X'])
                # Append all options that appear in the configuration file to the main command.
                for option in cl_options:
                    if isinstance(option, dict):
                        opt, val = option.items()[0]
                        cl.extend(['--{}'.format(opt), str(val)])
                    else:
                        cl.append('--{}'.format(option))

            #now add the base_mask for each lane
            for lane in sorted(lanes):
                #iterate thorugh each lane and add the correct --use-bases-mask for that lane
                #there is a single basemask for each lane, I checked it a couple of lines above
                base_mask = [per_lane_base_masks[lane][bm]['base_mask'] for bm in per_lane_base_masks[lane]][0] # get the base_mask
                base_mask_expr = "{}:".format(lane) + ",".join(base_mask)
                cl.extend(["--use-bases-mask", base_mask_expr])
        return cl

def _generate_clean_samplesheet(ssparser, fields_to_remove=None, rename_samples=True, rename_qPCR_suffix = False, fields_qPCR= None):
    """
        Will generate a 'clean' samplesheet, the given fields will be removed.
        if rename_samples is True, samples prepended with 'Sample_'  are renamed to match the sample name
        Will also replace 10X idicies like SI-GA-A3 with proper indicies like TGTGCGGG
    """
    output=""
    ##expand the ssparser if there are 10X lanes
    index_dict=parse_10X_indexes() #read the 10X indices
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

def look_for_lanes_with_10X_indicies(ssparser):
    """
    Given an ssparser object
    returns a list of lanes with 10X indicies
    """
    index_dict=parse_10X_indexes()
    tenX_lanes = set() #Set to only get each lane once
    not_tenX_lanes = set()
    for sample in ssparser.data:
        if sample['index'] in index_dict.keys():
            tenX_lanes.add(sample['Lane'])
        else:
            not_tenX_lanes.add(sample['Lane'])

    return (list(tenX_lanes),list(not_tenX_lanes))


def parse_10X_indexes():
    """
    Takes a file of 10X indexes and returns them as a dict.
    Todo: Set it up to take the file from config instead
    """
    index_dict={}
    indexfile="Chromium_10X_indexes.txt"
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
