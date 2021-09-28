

class Nanopore(object):
    """General Nanopore run"""
    def __init__(self):
        pass

    def is_not_transferred(run_id, transfer_log):
        """Return True if run id not in transfer.tsv, else False."""
        with open(transfer_log, 'r') as f:
            return run_id not in f.read()

    def transfer_run(run_dir):
        """rsync dir to Irma."""
        logger.info('Transferring run {} to analysis cluster'.format(run_dir))
        destination = CONFIG.get('nanopore_analysis').get('transfer').get('destination')
        rsync_opts = {'-Lav': None,
                    '--chown': ':ngi2016003',
                    '--chmod' : 'Dg+s,g+rw',
                    '-r' : None,
                    '--exclude' : 'work'}
        connection_details = CONFIG.get('nanopore_analysis').get('transfer').get('analysis_server')
        transfer_object = RsyncAgent(run_dir,
                                    dest_path=destination,
                                    remote_host=connection_details['host'],
                                    remote_user=connection_details['user'],
                                    validate=False,
                                    opts=rsync_opts)
        try:
            transfer_object.transfer()
        except RsyncError:
            logger.warn('An error occurred while transferring {} to the '
                        'ananlysis server. Please check the logfiles'.format(run_dir))
            return False
        return True

    def update_transfer_log(run_id, transfer_log):
        """Update transfer log with run id and date."""
        try:
            with open(transfer_log, 'a') as f:
                tsv_writer = csv.writer(f, delimiter='\t')
                tsv_writer.writerow([run_id, str(datetime.now())])
        except IOError:
            logger.warn('Could not update the transfer logfile for run {}. '
                        'Please make sure gets updated.'.format(run_id, transfer_log))
        return

    def archive_run(run_dir):
        """Move directory to nosync."""
        logger.info('Archiving run ' + run_dir)
        archive_dir = CONFIG.get('nanopore_analysis').get('finished_dir')
        top_dir = '/'.join(run_dir.split('/')[0:-2]) # Get the project folder to archive
        try:                                         # Try pathlib (pathlib.Path(run_dir).parent.parent) when running completely on python3
            shutil.move(top_dir, archive_dir)
            logger.info('Successfully archived {}'.format(run_dir))
        except shutil.Error:
            logger.warn('An error occurred when archiving {}. '
                        'Please check the logfile for more info.'.format(run_dir))
        return
