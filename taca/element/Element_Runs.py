import logging
import os

logger = logging.getLogger(__name__)


class Run:
    """Defines an Element run"""

    def __init__(self, run_dir, configuration):
        if not os.path.exists(run_dir):
            raise RuntimeError(f"Could not locate run directory {run_dir}")
        self.run_dir = os.path.abspath(run_dir)
        self.CONFIG = configuration
        self.demux_dir = "Demultiplexing"
    
    def is_transferred(self, transfer_file):
        pass
    
    def parse_rundir(self):
        pass