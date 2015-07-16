import os
import re
import csv
import glob
import datetime
import platform
import logging
from datetime import datetime
from flowcell_parser.classes import RunParametersParser, SampleSheetParser, RunParser

logger = logging.getLogger(__name__)

class Run(object):
    """ Defines an Illumina run
    """

    def __init__(self, run_dir, configuration):
        if not os.path.exists(run_dir) or not \
                os.path.exists(os.path.join(run_dir, 'runParameters.xml')):
            raise RuntimeError('Could not locate run directory {}'.format(run_dir))
        if 'analysis_server' not in configuration or \
            'bcl2fastq' not in configuration or \
                'samplesheets_dir' not in configuration:
            raise RuntimeError('configuration missing required entries (analysis_server, bcl2fastq, samplesheets_dir')
        self.run_dir = os.path.abspath(run_dir)
        self.id = os.path.basename(os.path.normpath(run_dir))
        pattern = r'(\d{6})_([ST-]*\w+\d+)_\d+_([AB]?)([A-Z0-9\-]+)'
        m = re.match(pattern, self.id)
        self.date        = m.group(1)
        self.instrument  = m.group(2)
        self.position    = m.group(3)
        self.flowcell_id = m.group(4)
        self.CONFIG      = configuration
        self._set_demux_folder(configuration)
        self.runParserObj = RunParser(self.run_dir)


    def demultiplex_run():
        raise NotImplementedError("Please Implement this method")

    def _set_sequencer_type(self, configuration):
        raise NotImplementedError("Please Implement this method")
    
    
    def _get_sequencer_type(self):
        if self.sequencer_type:
            return self.sequencer_type
        else:
            raise RuntimeError("sequencer_type not yet available!!")


    def _set_demux_folder(self, configuration):
        self.demux_dir = "Demultiplexing"
        for option in self.CONFIG['bcl2fastq']['options']:
            if isinstance(option, dict) and option.get('output-dir'):
                _demux_dir = option.get('output-dir')


    def _get_demux_folder(self):
        if self.demux_dir:
            return self.demux_dir
        else:
            raise RuntimeError("demux_folder not yet available!!")



    def _get_samplesheet(self):
        raise NotImplementedError("Please Implement this method")



    def _is_demultiplexing_done(self):
        if os.path.exists(os.path.join(self.run_dir, self._get_demux_folder(), 'Stats', 'DemultiplexingStats.xml')):
            return True
        else:
            return False



    def _is_demultiplexing_started(self):
        if os.path.exists(os.path.join(self.run_dir, self._get_demux_folder())):
            return True
        else:
            return False


    def _is_sequencing_done(self):
        if os.path.exists(os.path.join(self.run_dir, 'RTAComplete.txt')):
            return True
        else:
            return False


    def get_run_status(self):
        """
        return the status of the run, that is the trello card where it needs to be placed
        """
        demux_started   = self._is_demultiplexing_started() # True if demux is ongoing
        demux_done      = self._is_demultiplexing_done()    # True if demux is done
        sequencing_done = self._is_sequencing_done()        # True if sequencing is done
        if sequencing_done and demux_done:
            return 'COMPLETED' # run is done, tranfer is ongoing. Onother software/operator is responisble to move the run to nosync
        elif sequencing_done and demux_started and not demux_done:
            return 'IN_PROGRESS'
        elif sequencing_done and not demux_started:
            return 'TO_START'
        elif not sequencing_done:
            return 'SEQUENCING'
        else:
            raise RuntimeError('Unexpected status in get_run_status')




    def _generate_per_lane_base_mask(self):
        """
        This functions generate the base masks for each lane. This is required as often we
        run FC with lanes that have different index lengths. A more extreme case are lanes where
        multiple index lengths are present in the same lane.
        Hypotesis:
            - RunInfo.xml contains the configuration
            - this object contains a properly parsed samplesheet
            
        It returns an dict with a key for each lane:
        {lane1:
            {base_mask_string (e.g., Y150I6N2N8Y150):
                [ base_mask , [SampleSheetEntries]]
            }
            {base_mask_string (e.g., Y150I6N2I6N2Y150):
                [ base_mask , [SampleSheetEntries]]
            }
         lane2:
        }
            
        """
        #generate new ssparser (from the renamed smaplesheet)
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
            is_dual_index = False
            if data_entry.get('index'):
                index = data_entry['index']
                is_dual_index = "-" in  data_entry['index']
            elif data_entry.get('Index'):
                #specific for HiSeq, will disapper once we will use bcl2fastq_2.17
                index = data_entry['Index'].replace('-', '').replace('NoIndex', '')
                is_dual_index = "-" in  data_entry['Index']
            index_size = len(index)
            ##I have computed the index size, and I know if I have a single or dual index
            if is_dual_index:
                index_size = index_size/2
            #compute the basemask
            base_mask = self._compute_base_mask(runSetup, index_size, is_dual_index)
            base_mask_string = "".join(base_mask)
            #prepare the dictionary
            if base_mask_string not in base_masks[lane]:
                #first time I find such base mask in this lane,
                base_masks[lane][base_mask_string] = {'base_mask':base_mask,
                                                      'data' : []}
            base_masks[lane][base_mask_string]['data'].append(data_entry)

        return base_masks




    def _compute_base_mask(self, runSetup, index_size, dual_index_sample):
        """
            Assumptions:
                - if runSetup is of size 3, then single index run
                - if runSetup is of size 4, then dual index run
        """
        bm = []
        dual_index_run = False
        if len(runSetup) != 3 and len(runSetup) != 4:
            raise RuntimeError("when generating base_masks flooks like there are more than 3 or 4 reads in the RunSetup.xml")

        for read in runSetup:
            cycles = int(read['NumCycles'])
            if read['IsIndexedRead'] == 'N':
                bm.append('Y' + str(cycles))
            else:
                if index_size > cycles:
                    #the size of the index of the sample sheet is larger than the one specified by RunInfo.xml, somethig must be wrong
                    raise RuntimeError("when generating base_masks found index in samplesheet larger than the index specifed in RunInfo.xml")
                is_first_index_read = int(read['Number']) == 2
                #now prepare the base mask for this index read
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
                #when working on the second read index I need to know if the sample is dual index or not
                    if dual_index_sample:
                        i_remainder = cycles - index_size
                        if i_remainder > 0:
                            bm.append('I' + str(index_size) + 'N' + str(i_remainder))
                        else:
                            bm.append('I' + str(cycles))
                    else:
                    #if this sample is not dual index but the run is, then I need to ignore the second index completely
                            bm.append('N' + str(cycles))
        return bm
















