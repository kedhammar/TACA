import logging
import os
import json
from datetime import datetime

from taca.utils import misc
from taca.utils.filesystem import chdir

logger = logging.getLogger(__name__)


class Run:
    """Defines an Element run"""

    def __init__(self, run_dir, configuration):
        if not os.path.exists(run_dir):
            raise RuntimeError(f"Could not locate run directory {run_dir}")
        self.run_dir = os.path.abspath(run_dir)
        self.CONFIG = configuration
        self.demux_dir = os.path.join(self.run_dir, "Demultiplexing")
        self.final_sequencing_file = os.path.join(self.run_dir, "RunUploaded.json")
        self.demux_stats_file = os.path.join(self.demux_dir, "RunStats.json") #TODO: How to handle SideA/SideB?
        self.run_manifest_file = os.path.join(self.run_dir, "RunManifest.csv")
    
    def check_sequencing_status(self):
        if os.path.exists(self.final_sequencing_file):
            with open(self.final_sequencing_file) as json_file:
                sequencing_outcome = json.load(json_file).get("outcome")
            if sequencing_outcome != "OutcomeCompleted":
                return False
            else:
                return True
        else:
            return False
    
    def get_demultiplexing_status(self):
        if not os.path.exists(self.demux_dir):
            return "not started"
        elif os.path.exists(self.demux_dir) and not os.path.isfile(self.demux_stats_file):
            return "ongoing"
        elif os.path.exists(self.demux_dir) and os.path.isfile(self.demux_stats_file):
            return "finished"
    
    def manifest_exists(self):
        return os.path.isfile(self.run_manifest_file)
    
    def get_sample_info(self):
        sample_info = {} #TODO: populate 
        return sample_info
    
    def get_sample_types(self, sample_info):
        sample_types = () #TODO: populate 
        return sample_types
    
    def make_manifest(self, sample_info, sample_type):
        #TODO: make a manifest for a sample_type based on sample_info
        return
     
    def generate_demux_command(self):
        command = [self.CONFIG.get(self.software)["bin"], #TODO add path to bases2fastq executable to config
                   self.run_dir,
                   self.demux_dir, #TODO: how to handle SideA/SideB?
                   "-p 12"
                   ] 
        return command
    
    def start_demux(self):
        with chdir(self.run_dir):
            cmd = self.generate_demux_command()
            misc.call_external_command_detached(
                cmd, with_log_files=True, prefix=f"demux_"
            )
            logger.info(
                "Bases2Fastq conversion and demultiplexing "
                f"started for run {os.path.basename(self.id)} on {datetime.now()}"
            )
    
    def is_transferred(self, transfer_file):
        #TODO: return true if run in transfer log, else false
        pass
    
    def transfer_ongoing(self):
        #TODO: return true if hidden transfer file marker exists, else false
        pass
    
    def sync_metadata(self):
        #TODO: copy metadata from demuxed run to ngi-nas-ns
        pass
    
    def make_transfer_indicator(self):
        #TODO: touch a hidden file in the run directory
        pass
    
    def transfer(self):
        #TODO: rsync run to analysis cluster
        pass
    
    def remove_transfer_indicator(self):
        #TODO: remove hidden file in run directory
        pass
    
    def update_transfer_log(self, transfer_file):
        #TODO: update the transfer log
        pass
    
    def archive(self):
        #TODO: move run dir to nosync
        pass
    
    def parse_rundir(self):
        pass