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
            raise RunTimeError("sequencer_type not yet available!!")


    def _set_demux_folder(self, configuration):
        self.demux_dir = "Demultiplexing"
        for option in self.CONFIG['bcl2fastq']['options']:
            if isinstance(option, dict) and option.get('output-dir'):
                _demux_dir = option.get('output-dir')


    def _get_demux_folder(self):
        if self.demux_dir:
            return self.demux_dir
        else:
            raise RunTimeError("demux_folder not yet available!!")



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








