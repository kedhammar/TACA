import logging
import shutil

from taca.nanopore.nanopore import Nanopore
from taca.utils.config import CONFIG

logger = logging.getLogger(__name__)


class PromethION(Nanopore):
    """PromethION run"""
    def __init__(self, run_dir):
        super(PromethION, self).__init__(run_dir)
        self.transfer_log = CONFIG.get('nanopore_analysis').get('promethion_run').get('transfer').get('transfer_file')
        self.archive_dir = CONFIG.get('nanopore_analysis').get('promethion_run').get('finished_dir')
    
    def archive_run(self):
        """Move promethion directory to nosync."""
        logger.info('Archiving run ' + self.run_id)
        try:
            shutil.move(self.run_dir, self.archive_dir)
            logger.info('Successfully archived {}'.format(self.run_id))
            return True
        except shutil.Error:
            logger.warn('An error occurred when archiving {}. '
                        'Please check the logfile for more info.'.format(self.run_dir))
            return False
