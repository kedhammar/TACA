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
        self.summary_file = glob.glob(run_dir + '/final_summary*.txt')

    def is_not_transferred(self):
        """Return True if run id not in transfer.tsv, else False."""
        with open(self.transfer_log, 'r') as f:
            return self.run_id not in f.read()

    def transfer_run(self):
        """rsync dir to destination specified in config file."""
        destination = self.transfer_details.get('destination')
        rsync_opts = self.transfer_details.get('rsync_options')
        for k, v in rsync_opts.items():
            if v == 'None':
                rsync_opts[k] = None
        connection_details = self.transfer_details.get('analysis_server', None)
        logger.info('Transferring run {} to {}'.format(self.run_id, 
                                                       connection_details['host'] if connection_details 
                                                       else destination))
        if connection_details:
            transfer_object = RsyncAgent(self.run_dir,
                                        dest_path=destination,
                                        remote_host=connection_details.get('host'),
                                        remote_user=connection_details.get('user'),
                                        validate=False,
                                        opts=rsync_opts)
        else:
            transfer_object = RsyncAgent(self.run_dir,
                                        dest_path=destination,
                                        validate=False,
                                        opts=rsync_opts)
        try:
            transfer_object.transfer()
        except RsyncError:
            logger.warn('An error occurred while transferring {} to the '
                        'analysis server. Please check the logfiles'.format(self.run_dir))
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
        dir_to_move = str(pathlib.Path(self.run_dir).parent)
        top_dir = str(pathlib.Path(self.run_dir).parent.parent)  # Get the project folder to archive
        project_id = os.path.basename(top_dir)
        project_archive = os.path.join(self.archive_dir, project_id)
        if not os.path.exists(project_archive):
            os.mkdir(project_archive)
        try:
            print(dir_to_move, project_archive)
            shutil.move(dir_to_move, project_archive)
            logger.info('Successfully archived {}'.format(self.run_id))
            if not os.listdir(top_dir):
                logger.info("Project folder {} is empty. Removing it.".format(top_dir))
                os.rmdir(top_dir)
            else:
                logger.info("Some data is still left in {}. Keeping it.".format(top_dir))  # Might be another run for the same project
            return True
        except shutil.Error as e:
            logger.warn('The following error occurred when archiving {}:\n'
                        '{}'.format(self.run_dir, e))
            return False
