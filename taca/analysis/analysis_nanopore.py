"""Nanopore analysis methods for TACA."""
import os
import logging
import glob
import json
import html
import re

from dateutil.parser import parse
from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
from taca.utils.statusdb import NanoporeRunsConnection
from taca.nanopore.minion import MinIONdelivery, MinIONqc
from taca.nanopore.ont_transfer import PromethionTransfer, MinionTransfer

logger = logging.getLogger(__name__)



def is_date(string):
    """
    Return whether the string can be interpreted as a date.

    :param string: str, string to check for date
    From https://stackoverflow.com/questions/25341945/check-if-string-has-date-any-format
    """
    try: 
        parse(string, fuzzy=False)
        return True
    except ValueError:
        return False

def find_minion_runs(minion_data_dir, skip_dirs):
    """Find nanopore runs to process."""
    found_run_dirs = []
    try:
        found_top_dirs = [os.path.join(minion_data_dir, top_dir) for top_dir in os.listdir(minion_data_dir)
                            if os.path.isdir(os.path.join(minion_data_dir, top_dir))
                            and top_dir not in skip_dirs]
    except OSError:
        logger.warning('There was an issue locating the following directory: {}. '
                    'Please check that it exists and try again.'.format(minion_data_dir))
    # Get the actual location of the run directories in /var/lib/MinKnow/data/QC_runs/USERDETERMINEDNAME/USERDETSAMPLENAME/run
    if found_top_dirs:
        for top_dir in found_top_dirs:
            if os.path.isdir(top_dir):
                for sample_dir in os.listdir(top_dir):
                    if os.path.isdir(os.path.join(top_dir, sample_dir)):
                        for run_dir in os.listdir(os.path.join(top_dir, sample_dir)):
                            found_run_dirs.append(os.path.join(top_dir, sample_dir, run_dir))
    else:
        logger.warning('Could not find any run directories in {}'.format(minion_data_dir))
    return found_run_dirs

def find_ont_transfer_runs(ont_data_dir, skip_dirs):
    """Find runs in ngi-nas. 
    These are assumed to be flowcell dirs, not project dirs.
    """
    try:
        found_dirs = [os.path.join(ont_data_dir, top_dir) for top_dir in os.listdir(ont_data_dir)
                            if os.path.isdir(os.path.join(ont_data_dir, top_dir))
                            and top_dir not in skip_dirs]
    except OSError:
        logger.warning('There was an issue locating the following directory: {}. '
                    'Please check that it exists and try again.'.format(ont_data_dir))
    return found_dirs
        

def check_ongoing_sequencing(run_dir):
    """Check if sequencing is ongoing for the given run dir """
    summary_file = glob.glob(run_dir + '/final_summary*.txt')

    if not len(summary_file):
        logger.info('Sequencing ongoing for run {}. Will not start nanoseq at this time'.format(run_dir))
        return True
    else:
        return False

def process_minion_qc_run(minion_run, sequencing_ongoing=False, nanoseq_ongoing=False):
    """Process MinION QC runs.
    
    Will not start nanoseq if a sequencing run is ongoing, to limit memory usage.
    Will also maximum start one nanoseq run at once, for the same reason.
    """
    logger.info('Processing QC run: {}'.format(minion_run.run_dir))
    email_recipients = CONFIG.get('mail').get('recipients')

    if len(minion_run.summary_file) and os.path.isfile(minion_run.summary_file[0]) and not os.path.isdir(minion_run.nanoseq_dir):
        logger.info('Sequencing done for run {}. Attempting to start analysis.'.format(minion_run.run_dir))
        if not minion_run.nanoseq_sample_sheet:
            minion_run.parse_lims_sample_sheet()

        if os.path.isfile(minion_run.nanoseq_sample_sheet):
            if nanoseq_ongoing:
                logger.warning('Nanoseq already started, will not attempt to start for {}'.format(minion_run.run_dir))
            elif sequencing_ongoing:
                logger.warning('Sequencing ongoing, will not attempt to start Nanoseq for {}'.format(minion_run.run_dir))
            else:
                minion_run.start_nanoseq()
                nanoseq_ongoing = True

        else:
            logger.warning('Samplesheet not found for run {}. Operator notified. Skipping.'.format(minion_run.run_dir))
            email_subject = ('Samplesheet missing for run {}'.format(os.path.basename(minion_run.run_dir)))
            email_message = 'There was an issue locating the samplesheet for run {}.'.format(minion_run.run_dir)
            send_mail(email_subject, email_message, email_recipients)

    elif os.path.isdir(minion_run.nanoseq_dir) and not os.path.isfile(minion_run.nanoseq_exit_status_file):
        logger.info('Nanoseq has started for run {} but is not yet done. Skipping.'.format(minion_run.run_dir))

    elif os.path.isdir(minion_run.nanoseq_dir) and os.path.isfile(minion_run.nanoseq_exit_status_file):
        nanoseq_successful = minion_run.check_exit_status(minion_run.nanoseq_exit_status_file)
        if nanoseq_successful:
            if not os.path.isdir(minion_run.anglerfish_dir):
                logger.info('Nanoseq done for run {}. Attempting to start Anglerfish.'.format(minion_run.run_id))
                if not minion_run.anglerfish_sample_sheet:
                    minion_run.anglerfish_sample_sheet = os.path.join(minion_run.run_dir, 'anglerfish_sample_sheet.csv')  # For cronjob, AF sample sheet was generated at previous run
                if os.path.isfile(minion_run.anglerfish_sample_sheet):
                    minion_run.start_anglerfish()
                else:
                    logger.warning('Anglerfish sample sheet missing for run {}. '
                                'Please provide one using --anglerfish_sample_sheet '
                                'if running TACA manually.'.format(minion_run.run_id))

            elif not os.path.isfile(minion_run.anglerfish_exit_status_file):
                logger.info('Anglerfish has started for run {} but is not yet done. Skipping.'.format(minion_run.run_id))

            elif os.path.isfile(minion_run.anglerfish_exit_status_file):
                anglerfish_successful = minion_run.check_exit_status(minion_run.anglerfish_exit_status_file)
                if anglerfish_successful:
                    if minion_run.copy_results_for_lims():
                        logger.info('Anglerfish finished OK for run {}. Notifying operator.'.format(minion_run.run_id))
                        email_subject = ('Anglerfish successfully processed run {}'.format(minion_run.run_id))
                        email_message = ('Anglerfish has successfully finished for run {}. Please '
                                         'finish the QC step in lims.').format(minion_run.run_id)
                        send_mail(email_subject, email_message, email_recipients)
                    else:
                        email_subject = ('Run processed with errors: {}'.format(minion_run.run_id))
                        email_message = ('Anglerfish has successfully finished for run {} but an error '
                                         'occurred while transferring the results to lims.').format(minion_run.run_id)
                        send_mail(email_subject, email_message, email_recipients)

                    if minion_run.is_not_transferred():
                        if minion_run.transfer_run():
                            if minion_run.update_transfer_log():
                                logger.info('Run {} has been synced to the analysis cluster.'.format(minion_run.run_id))
                            else:
                                email_subject = ('Run processed with errors: {}'.format(minion_run.run_id))
                                email_message = ('Run {} has been transferred, but an error occurred while updating '
                                                 'the transfer log').format(minion_run.run_id)
                                send_mail(email_subject, email_message, email_recipients)

                            if minion_run.archive_run():
                                logger.info('Run {} is finished and has been archived. Notifying operator.'.format(minion_run.run_id))
                                email_subject = ('Run successfully processed: {}'.format(minion_run.run_id))
                                email_message = ('Run {} has been analysed, transferred and archived '
                                                 'successfully.').format(minion_run.run_id)
                                send_mail(email_subject, email_message, email_recipients)
                            else:
                                email_subject = ('Run processed with errors: {}'.format(minion_run.run_id))
                                email_message = ('Run {} has been analysed, but an error occurred during '
                                                 'archiving').format(minion_run.run_id)
                                send_mail(email_subject, email_message, email_recipients)

                        else:
                            logger.warning('An error occurred during transfer of run {} '
                                        'to the analysis cluster. Notifying operator.'.format(minion_run.run_id))
                            email_subject = ('Run processed with errors: {}'.format(minion_run.run_id))
                            email_message = ('Run {} has been analysed, but an error occurred during '
                                             'transfer.').format(minion_run.run_id)
                            send_mail(email_subject, email_message, email_recipients)

                    else:
                        logger.warning('The following run has already been transferred, '
                                    'skipping: {}'.format(minion_run.run_id))

                else:
                    logger.warning('Anglerfish exited with a non-zero exit status for run {}. '
                                'Notifying operator.'.format(minion_run.run_id))
                    email_subject = ('Run processed with errors: {}'.format(minion_run.run_id))
                    email_message = ('Anglerfish exited with errors for run {}. Please '
                                     'check the log files and restart.').format(minion_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)

        else:
            logger.warning('Nanoseq exited with a non-zero exit status for run {}. '
                        'Notifying operator.'.format(minion_run.run_id))
            email_subject = ('Analysis failed for run {}'.format(minion_run.run_id))
            email_message = 'The nanoseq analysis failed for run {}.'.format(minion_run.run_id)
            send_mail(email_subject, email_message, email_recipients)

    else:
        logger.info('Run {} not finished sequencing yet. Skipping.'.format(minion_run.run_id))

    return nanoseq_ongoing

def process_minion_delivery_run(minion_run):
    """Process minion delivery runs."""
    email_recipients = CONFIG.get('mail').get('recipients')
    logger.info('Processing run {}'.format(minion_run.run_id))
    minion_run.dump_path()
    if not len(minion_run.summary_file):  # Run not finished, only rsync
        minion_run.transfer_run()
    else:  # Run finished, rsync and archive
        if minion_run.transfer_run():
            minion_run.archive_run()
            logger.info('Run {} has been fully transferred.'.format(minion_run.run_id))
            email_subject = ('Run successfully processed: {}'.format(minion_run.run_id))
            email_message = ('Run {} has been transferred and archived '
                            'successfully.').format(minion_run.run_id)
            send_mail(email_subject, email_message, email_recipients)
        else:
            logger.warning('An error occurred during transfer of run {}.'.format(minion_run.run_id))
            email_subject = ('Run processed with errors: {}'.format(minion_run.run_id))
            email_message = ('An error occurred during the '
                            'transfer of run {}.').format(minion_run.run_id)
            send_mail(email_subject, email_message, email_recipients)
            

def ont2couch(ont_run):
    """ Check run vs statusdb. Create or update run entry as needed.
    """

    run_pattern = re.compile("^(\d{8})_(\d{4})_([0-9a-zA-Z]+)_([0-9a-zA-Z]+)_([0-9a-zA-Z]+)$")

    sesh = NanoporeRunsConnection(CONFIG["statusdb"], dbname="nanopore_runs")        

    if re.match(run_pattern, ont_run.run_id):
        logger.debug(f"Run {ont_run.run_id} looks like a run directory, continuing.")
    else:
        error_message = f"Run {ont_run.run_id} does not match the regex of a run directory (yyyymmdd_hhmm_pos|device_fcID_hash)."
        logger.error(error_message)
        raise AssertionError(error_message)
    
    # If no run document exists in the database
    if not sesh.check_run_exists(ont_run):
        logger.debug(f"Run {ont_run.run_id} does not exist in the database, creating entry for ongoing run.")
        
        # Create an ongoing run document
        run_path_file = os.path.join(ont_run.run_dir, 'run_path.txt')
        sesh.create_ongoing_run(ont_run, open(run_path_file, "r").read().strip())
        logger.debug(f"Successfully created db entry for ongoing run {ont_run.run_id}.")
        
    # If a run document DOES exist in the database
    else:
        # If the run document is marked as "ongoing"
        if sesh.check_run_status(ont_run) == "ongoing":
            
            logger.debug(f"Run {ont_run.run_id} exists in the database as an ongoing run.")

            # If the run is finished
            if len(ont_run.summary_file) != 0:

                logger.debug(f"Run {ont_run.run_id} has finished sequencing, updating the db entry.")

                # Parse the MinKNOW .json report file and finish the ongoing run document
                glob_json = glob.glob(ont_run.run_dir + '/report*.json')
                glob_html = glob.glob(ont_run.run_dir + '/report*.html')

                if len(glob_json) == 0 or len(glob_html) == 0:
                    error_message = f"Run {ont_run.run_id} is marked as finished, but missing report files."
                    logger.error(error_message)
                    raise AssertionError(error_message)
                
                dict_json = json.load(open(glob_json[0], "r"))
                dict_html = {
                    "minknow_report_name": glob_html[0].split("/")[-1],
                    "minknow_report_content": html.escape(open(glob_html[0], "r").read())
                }

                sesh.finish_ongoing_run(ont_run, dict_json, dict_html)
                logger.debug(f"Successfully updated the db entry of run {ont_run.run_id}")

                # Transfer the MinKNOW .html report file, may require SSH access via IT
                scp_command = f'scp {glob_html[0]} {os.path.join(CONFIG["nanopore_analysis"]["ont_transfer"]["minknow_reports_dir"], glob_html[0])}'
                os.system(scp_command)
                logger.debug(f"Successfully transferred the MinKNOW report of run {ont_run.run_id}")

            else:
                logger.debug(f"Run {ont_run.run_id} has not finished sequencing, do nothing.")
        
        else:
            logger.debug(f"Run {ont_run.run_id} exists in the database as an finished run, do nothing.")


def transfer_ont_run(ont_run):
    """Transfer ONT runs to HPC cluster."""
    email_recipients = CONFIG.get('mail').get('recipients')
    logger.info('Processing run {}'.format(ont_run.run_id))

    # Update StatusDB
    try:
        ont2couch(ont_run)
        logger.info(f"Database update for run {ont_run.run_id} successful")
    except Exception as e:
        logger.warning(f"Database update for run {ont_run.run_id} failed")
        email_subject = ('Run processed with errors: {}'.format(ont_run.run_id))
        email_message = (f"An error occured when updating statusdb with run {ont_run.run_id}.\n{e}")
        send_mail(email_subject, email_message, email_recipients)

    if len(ont_run.summary_file) and os.path.isfile(ont_run.summary_file[0]):
        logger.info('Sequencing done for run {}. Attempting to start processing.'.format(ont_run.run_id))
        if ont_run.is_not_transferred():
            if ont_run.transfer_run():
                if ont_run.update_transfer_log():
                    logger.info('Run {} has been synced to the analysis cluster.'.format(ont_run.run_id))
                else:
                    email_subject = ('Run processed with errors: {}'.format(ont_run.run_id))
                    email_message = ('Run {} has been transferred, but an error occurred while updating '
                                        'the transfer log').format(ont_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)
            
                if ont_run.archive_run():
                    logger.info('Run {} is finished and has been archived. '
                                'Notifying operator.'.format(ont_run.run_id))
                    email_subject = ('Run successfully processed: {}'.format(ont_run.run_id))
                    email_message = ('Run {} has been transferred and archived '
                                    'successfully.').format(ont_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)
                else:
                    email_subject = ('Run processed with errors: {}'.format(ont_run.run_id))
                    email_message = ('Run {} has been analysed, but an error occurred during '
                                    'archiving').format(ont_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)

            else:
                email_subject = ('Run processed with errors: {}'.format(ont_run.run_id))
                email_message = ('An error occurred during transfer of run {} '
                                'to the analysis cluster.').format(ont_run.run_id)
                send_mail(email_subject, email_message, email_recipients)

        else:
            logger.warning('The following run has already been transferred, '
                        'skipping: {}'.format(ont_run.run_id))
    else:
        logger.info('Run {} not finished sequencing yet. Skipping.'.format(ont_run.run_id))


def process_minion_qc_runs(run, nanoseq_sample_sheet, anglerfish_sample_sheet):
    """Find and process MinION QC runs on Squiggle."""
    if run:
        if is_date(os.path.basename(run).split('_')[0]):
            minion_run = MinIONqc(os.path.abspath(run), nanoseq_sample_sheet, anglerfish_sample_sheet)
            process_minion_qc_run(minion_run)
        else:
            logger.warning('The specified path is not a flow cell. Please '
                        'provide the full path to the flow cell you wish to process.')
    else:
        nanopore_data_dir = CONFIG.get('nanopore_analysis').get('minion_qc_run').get('data_dir')
        skip_dirs = CONFIG.get('nanopore_analysis').get('minion_qc_run').get('ignore_dirs')
        runs_to_process = find_minion_runs(nanopore_data_dir, skip_dirs)
        sequencing_ongoing, nanoseq_ongoing = False, False
        for run_dir in runs_to_process:
            if check_ongoing_sequencing(run_dir):
                sequencing_ongoing = True
                break
        for run_dir in runs_to_process:
            minion_run = MinIONqc(run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet)
            nanoseq_ongoing = process_minion_qc_run(minion_run,
                                                    sequencing_ongoing=sequencing_ongoing,
                                                    nanoseq_ongoing=nanoseq_ongoing)

def process_minion_delivery_runs(run):
    """Find MinION delivery runs on Squiggle and transfer them to ngi-nas."""
    if run:
        if is_date(os.path.basename(run).split('_')[0]):
            minion_run = MinIONdelivery(os.path.abspath(run))
            process_minion_delivery_run(minion_run)
        else:
            logger.warning('The specified path is not a flow cell. Please '
                        'provide the full path to the flow cell you wish to process.')
    else:
        minion_data_dir = CONFIG.get('nanopore_analysis').get('minion_delivery_run').get('data_dir')
        skip_dirs = CONFIG.get('nanopore_analysis').get('minion_delivery_run').get('ignore_dirs')
        runs_to_process = find_minion_runs(minion_data_dir, skip_dirs)
        for run_dir in runs_to_process:
            minion_run = MinIONdelivery(run_dir)
            process_minion_delivery_run(minion_run)

def transfer_finished(run):
    """Find finished ONT runs in ngi-nas and transfer to HPC cluster."""
    if run:
        if is_date(os.path.basename(run).split('_')[0]):
            if 'minion' in run:
                ont_run = MinionTransfer(os.path.abspath(run))
                transfer_ont_run(ont_run)
            elif 'promethion' in run:
                ont_run = PromethionTransfer(os.path.abspath(run))
                transfer_ont_run(ont_run)
        else:
            logger.warning('The specified path is not a flow cell. Please '
                        'provide the full path to the flow cell you wish to process.')
    else:
        # Locate all runs in /srv/ngi_data/sequencing/promethion and /srv/ngi_data/sequencing/minion
        ont_data_dirs = CONFIG.get('nanopore_analysis').get('ont_transfer').get('data_dirs')
        skip_dirs = CONFIG.get('nanopore_analysis').get('ont_transfer').get('ignore_dirs')
        for data_dir in ont_data_dirs:
            runs_to_process = find_ont_transfer_runs(data_dir, skip_dirs)
            for run_dir in runs_to_process:
                if 'minion' in data_dir:
                    ont_run = MinionTransfer(run_dir)
                    transfer_ont_run(ont_run)
                elif 'promethion' in data_dir:
                    ont_run = PromethionTransfer(run_dir)
                    transfer_ont_run(ont_run)
