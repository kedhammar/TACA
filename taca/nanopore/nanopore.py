import os
import logging
import csv
import shutil
import glob
import pathlib

from datetime import datetime
from taca.utils.config import CONFIG
from taca.utils.transfer import RsyncAgent, RsyncError

logger = logging.getLogger(__name__)

class Nanopore(object):
    """General Nanopore run"""
    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.run_id = os.path.basename(run_dir)
        self.transfer_log = CONFIG.get('nanopore_analysis').get('transfer').get('transfer_file')
        self.summary_file = glob.glob(run_dir + '/final_summary*.txt')

    def is_not_transferred(self):
        """Return True if run id not in transfer.tsv, else False."""
        with open(self.transfer_log, 'r') as f:
            return self.run_id not in f.read()

    def transfer_run(self):
        """rsync dir to Irma."""
        logger.info('Transferring run {} to analysis cluster'.format(self.run_id))
        destination = CONFIG.get('nanopore_analysis').get('transfer').get('destination')
        rsync_opts = CONFIG.get('nanopore_analysis').get('transfer').get('rsync_options')
        for k, v in rsync_opts.items():
            if v == 'None':
                rsync_opts[k] = None
        connection_details = CONFIG.get('nanopore_analysis').get('transfer').get('analysis_server')
        transfer_object = RsyncAgent(self.run_dir,
                                    dest_path=destination,
                                    remote_host=connection_details['host'],
                                    remote_user=connection_details['user'],
                                    validate=False,
                                    opts=rsync_opts)
        try:
            transfer_object.transfer()
        except RsyncError:
            logger.warn('An error occurred while transferring {} to the '
                        'ananlysis server. Please check the logfiles'.format(self.run_dir))
            return False
        return True

    def update_transfer_log(self):
        """Update transfer log with run id and date."""
        try:
            with open(self.transfer_log, 'a') as f:
                tsv_writer = csv.writer(f, delimiter='\t')
                tsv_writer.writerow([self.run_id, str(datetime.now())])
                return True
        except IOError:
            logger.warn('Could not update the transfer logfile for run {}. '
                        'Please make sure it gets updated.'.format(self.run_id, self.transfer_log))
            return False

    def archive_run(self):
        """Move directory to nosync."""
        logger.info('Archiving run ' + self.run_id)
        archive_dir = CONFIG.get('nanopore_analysis').get('finished_dir')
        top_dir = str(pathlib.Path(self.run_dir).parent.parent)  # Get the project folder to archive
        try:
            shutil.move(top_dir, archive_dir)
            logger.info('Successfully archived {}'.format(self.run_id))
            return True
        except shutil.Error:
            logger.warn('An error occurred when archiving {}. '
                        'Please check the logfile for more info.'.format(self.run_dir))
            return False
