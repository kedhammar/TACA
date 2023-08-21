import logging
import shutil
import os

from taca.nanopore.nanopore import Nanopore
from taca.utils.config import CONFIG

logger = logging.getLogger(__name__)


class ONTTransfer(Nanopore):
    """Base class for transfer of ONT data to HPC cluster"""
    def __init__(self, run_dir):
        super(ONTTransfer, self).__init__(run_dir)
        self.sync_finished_indicator = os.path.join(run_dir, '.sync_finished')
    
    def archive_run(self):
        """Move run directory to nosync."""
        logger.info('Archiving run ' + self.run_id)
        try:
            shutil.move(self.run_dir, self.archive_dir)
            logger.info('Successfully archived {}'.format(self.run_id))
            return True
        except shutil.Error:
            logger.warn('An error occurred when archiving {}. '
                        'Please check the logfile for more info.'.format(self.run_dir))
            return False

class PromethionTransfer(ONTTransfer):
    """Class for transfer of PromethION data to HPC cluster"""
    def __init__(self, run_dir):
        super(PromethionTransfer, self).__init__(run_dir)
        self.transfer_details = CONFIG.get('nanopore_analysis').get('ont_transfer').get('promethion')
        self.transfer_log = self.transfer_details.get('transfer_file')
        self.archive_dir = self.transfer_details.get('finished_dir')
        

class MinionTransfer(ONTTransfer):
    """Class for transfer of MinION data to HPC cluster"""
    def __init__(self, run_dir):
        super(MinionTransfer, self).__init__(run_dir)
        self.transfer_details = CONFIG.get('nanopore_analysis').get('ont_transfer').get('minion')
        self.transfer_log = self.transfer_details.get('transfer_file')
        self.archive_dir = self.transfer_details.get('finished_dir')