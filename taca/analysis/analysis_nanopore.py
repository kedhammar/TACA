"""Nanopore analysis methods for TACA."""
import os
import logging
import glob
import csv
import subprocess
import shutil
import smtplib

from datetime import datetime
from taca.utils.config import CONFIG
from taca.utils.minion_barcodes import BARCODES
from taca.utils.transfer import RsyncAgent, RsyncError
from taca.utils.misc import send_mail

logger = logging.getLogger(__name__)

def find_runs_to_process():
    """Find nanopore runs to process."""
    nanopore_data_dir = CONFIG.get('nanopore_analysis').get('data_dir')[0]
    found_run_dirs = []
    try:
        found_top_dirs = [os.path.join(nanopore_data_dir, top_dir) for top_dir in os.listdir(nanopore_data_dir)
                 if os.path.isdir(os.path.join(nanopore_data_dir, top_dir))
                 and top_dir != 'nosync']
    except OSError:
        logger.warn('There was an issue locating the following directory: {}. '
                    'Please check that it exists and try again.'.format(nanopore_data_dir))
    # Get the actual location of the run directories in /var/lib/MinKnow/data/USERDETERMINEDNAME/USERDETSAMPLENAME/run
    if found_top_dirs:
        for top_dir in found_top_dirs:
            for sample_dir in os.listdir(top_dir):
                for run_dir in os.listdir(os.path.join(top_dir, sample_dir)):
                    found_run_dirs.append(os.path.join(top_dir, sample_dir, run_dir))
    else:
        logger.warn('Could not find any run directories in {}'.format(nanopore_data_dir))
    return found_run_dirs

def process_minion_run(run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet):
    """Process minion runs."""
    qc_run = True
    if nanoseq_sample_sheet and not anglerfish_sample_sheet:
        qc_run = False

    logger.info('Processing run: {} as a {}'.format(run_dir, 'QC run' if qc_run else 'non-QC run'))
    summary_file = glob.glob(run_dir + '/final_summary*.txt')[0]
    nanoseq_dir = os.path.join(run_dir, 'nanoseq_output')
    anglerfish_dir = os.path.join(run_dir, 'anglerfish_output')
    anglerfish_sample_sheet = os.path.join(run_dir, 'anglerfish_sample_sheet.csv')
    nanoseq_exit_status_file = os.path.join(run_dir, '.exitcode_for_nanoseq')
    anglerfish_exit_status_file = os.path.join(run_dir, '.exitcode_for_anglerfish')
    email_recipients = CONFIG.get('mail').get('recipients')

    if os.path.isfile(summary_file) and not os.path.isdir(nanoseq_dir):
        logger.info('Sequencing done for run {}. Attempting to start analysis.'.format(run_dir))
        if not nanoseq_sample_sheet:
            nanoseq_sample_sheet = parse_lims_sample_sheet(run_dir)

        if os.path.isfile(nanoseq_sample_sheet):
            start_nanoseq(run_dir, nanoseq_sample_sheet)

        else:
            logger.warn('Samplesheet not found for run {}. Operator notified. Skipping.'.format(run_dir))
            email_subject = ('Samplesheet missing for run {}'.format(os.path.basename(run_dir)))
            email_message = 'There was an issue locating the samplesheet for run {}.'.format(run_dir)
            send_mail(email_subject, email_message, email_recipients)

    elif os.path.isdir(nanoseq_dir) and not os.path.isfile(nanoseq_exit_status_file):
        logger.info('Nanoseq has started for run {} but is not yet done. Skipping.'.format(run_dir))

    elif os.path.isdir(nanoseq_dir) and os.path.isfile(nanoseq_exit_status_file):
        nanoseq_successful = check_exit_status(nanoseq_exit_status_file)
        if nanoseq_successful:
            run_id = os.path.basename(run_dir)
            transfer_log = CONFIG.get('nanopore_analysis').get('transfer').get('transfer_file')

            if qc_run and not os.path.isdir(anglerfish_dir):
                logger.info('Nanoseq done for run {}. Attempting to start Anglerfish.'.format(run_id))
                start_anglerfish(run_dir, anglerfish_sample_sheet, anglerfish_dir)

            elif qc_run and not os.path.isfile(anglerfish_exit_status_file):
                logger.info('Anglerfish has started for run {} but is not yet done. Skipping.'.format(run_id))

            elif qc_run and os.path.isfile(anglerfish_exit_status_file):
                anglerfish_successful = check_exit_status(anglerfish_exit_status_file)
                if anglerfish_successful:
                    copy_results_for_lims(run_dir, anglerfish_dir)
                    logger.info('Anglerfish finished OK for run {}. Notifying operator.'.format(run_id))
                    email_subject = ('Anglerfish successfully processed run {}'.format(os.path.basename(run_id)))
                    email_message = ('Anglerfish has successfully finished for run {}. Please '
                                     'finish the QC step in lims.').format(run_id)
                    send_mail(email_subject, email_message, email_recipients)

                    if is_not_transferred(run_id, transfer_log):
                        if transfer_run(run_dir):
                            update_transfer_log(run_id, transfer_log)
                            logger.info('Run {} has been synced to the analysis cluster.'.format(run_id))
                            archive_run(run_dir)
                            logger.info('Run {} is finished and has been archived. Notifying operator.'.format(run_id))
                            email_subject = ('Run successfully processed: {}'.format(os.path.basename(run_id)))
                            email_message = ('Run {} has been analysed, transferred and archived '
                                             'successfully.').format(run_id)
                            send_mail(email_subject, email_message, email_recipients)

                        else:
                            logger.warn('An error occurred during transfer of run {} '
                                        'to Irma. Notifying operator.'.format(run_dir))
                            email_subject = ('Run processed with errors: {}'.format(os.path.basename(run_id)))
                            email_message = ('Run {} has been analysed, but an error occurred during '
                                             'transfer.').format(run_id)
                            send_mail(email_subject, email_message, email_recipients)

                    else:
                        logger.warn('The following run has already been transferred, '
                                    'skipping: {}'.format(run_dir))

                else:
                    logger.warn('Anglerfish exited with a non-zero exit status for run {}. '
                                'Notifying operator.'.format(run_dir))
                    email_subject = ('Run processed with errors: {}'.format(os.path.basename(run_id)))
                    email_message = ('Anglerfish exited with errors for run {}. Please '
                                     'check the log files and restart.').format(run_id)
                    send_mail(email_subject, email_message, email_recipients)

            elif not qc_run:
                if is_not_transferred(run_id, transfer_log):
                    if transfer_run(run_dir):
                        update_transfer_log(run_id, transfer_log)
                        logger.info('Run {} has been synced to the analysis cluster.'.format(run_id))
                        archive_run(run_dir)
                        logger.info('Run {} is finished and has been archived. Notifying operator.'.format(run_id))
                        email_subject = ('Run successfully processed: {}'.format(run_id))
                        email_message = ('Run {} has been analysed, transferred and archived '
                                         'successfully.').format(run_id)
                        send_mail(email_subject, email_message, email_recipients)

                    else:
                        logger.warn('An error occurred during transfer of run {} '
                                    'to Irma. Notifying operator.'.format(run_dir))
                        email_subject = ('Run processed with errors: {}'.format(run_id))
                        email_message = ('Run {} has been analysed, but an error occurred during '
                                         'transfer.').format(run_id)
                        send_mail(email_subject, email_message, email_recipients)

                else:
                    logger.warn('The following run has already been transferred, '
                                'skipping: {}'.format(run_id))

        else:
            logger.warn('Nanoseq exited with a non-zero exit status for run {}. '
                        'Notifying operator.'.format(run_dir))
            email_subject = ('Analysis failed for run {}'.format(os.path.basename(run_dir)))
            email_message = 'The nanoseq analysis failed for run {}.'.format(run_dir)
            send_mail(email_subject, email_message, email_recipients)

    else:
        logger.info('Run {} not finished sequencing yet. Skipping.'.format(run_dir))

    return

    def process_promethion_run(run_dir):
        """Process promethion runs."""
        pass


def process_minion_runs(run, nanoseq_sample_sheet, anglerfish_sample_sheet):
    """Find runs and kick off processing."""
    if run:
        process_minion_run(os.path.abspath(run), nanoseq_sample_sheet, anglerfish_sample_sheet)
    else:
        runs_to_process = find_runs_to_process()
        for run_dir in runs_to_process:
            process_minion_run(run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet)

def process_promethion_runs(run):
    """Find runs and kick off processing."""
    if run:
        process_promethion_run(os.path.abspath(run))
    else:
        runs_to_process = find_runs_to_process()
        for run_dir in runs_to_process:
            process_promethion_run(run_dir)