import logging
import shutil

from taca.nanopore.nanopore import Nanopore
from taca.utils.config import CONFIG

logger = logging.getLogger(__name__)


class ONTTransfer(Nanopore):
    """ONT run for transfer to HPC cluster"""
    def __init__(self, run_dir):
        super(ONTTransfer, self).__init__(run_dir)
    
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
    """PromethION run for transfer to HPC cluster"""
    def __init__(self, run_dir):
        super(PromethionTransfer, self).__init__(run_dir)
        self.transfer_details = CONFIG.get('nanopore_analysis').get('ont_transfer').get('promethion_transfer')
        self.transfer_log = self.transfer_details.get('transfer_file') #TODO: update config
        self.archive_dir = self.transfer_details.get('finished_dir')
        

class MinionTransfer(ONTTransfer):
    """MinION run for transfer to HPC cluster"""
    def __init__(self, run_dir):
        super(MinionTransfer, self).__init__(run_dir)
        self.transfer_details = CONFIG.get('nanopore_analysis').get('ont_transfer').get('minion_transfer')
        self.transfer_log = self.transfer_details.get('transfer_file') #TODO: update config
        self.archive_dir = self.transfer_details.get('finished_dir')