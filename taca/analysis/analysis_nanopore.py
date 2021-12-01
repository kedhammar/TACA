"""Nanopore analysis methods for TACA."""
import os
import logging
import glob

from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
from taca.nanopore.minion import MinION
from taca.nanopore.promethion import PromethION

logger = logging.getLogger(__name__)

def find_runs_to_process():
    """Find nanopore runs to process."""
    nanopore_data_dir = CONFIG.get('nanopore_analysis').get('data_dir')[0]
    found_run_dirs = []
    skip_dirs = CONFIG.get('nanopore_analysis').get('ignore_dirs')
    try:
        found_top_dirs = [os.path.join(nanopore_data_dir, top_dir) for top_dir in os.listdir(nanopore_data_dir)
                            if os.path.isdir(os.path.join(nanopore_data_dir, top_dir))
                            and top_dir not in skip_dirs]
    except OSError:
        logger.warn('There was an issue locating the following directory: {}. '
                    'Please check that it exists and try again.'.format(nanopore_data_dir))
    # Get the actual location of the run directories in /var/lib/MinKnow/data/QC_runs/USERDETERMINEDNAME/USERDETSAMPLENAME/run
    if found_top_dirs:
        for top_dir in found_top_dirs:
            if os.path.isdir(top_dir):
                for sample_dir in os.listdir(top_dir):
                    if os.path.isdir(os.path.join(top_dir, sample_dir)):
                        for run_dir in os.listdir(os.path.join(top_dir, sample_dir)):
                            found_run_dirs.append(os.path.join(top_dir, sample_dir, run_dir))
    else:
        logger.warn('Could not find any run directories in {}'.format(nanopore_data_dir))
    return found_run_dirs

def check_ongoing_sequencing(run_dir):
    """Check if sequencing is ongoing for the given run dir """
    summary_file = glob.glob(run_dir + '/final_summary*.txt')

    if not len(summary_file):
        logger.info('Sequencing ongoing for run {}. Will not start nanoseq at this time'.format(run_dir))
        return True
    else:
        return False

def process_minion_run(minion_run, sequencing_ongoing=False, nanoseq_ongoing=False):
    """Process MinION QC runs.
    
    Will not start nanoseq if a sequencing run is ongoing, to limit memory usage.
    Will also maximum start one nanoseq run at once, for the same reason.
    """
    logger.info('Processing run: {}'.format(minion_run.run_dir))
    email_recipients = CONFIG.get('mail').get('recipients')

    if len(minion_run.summary_file) and os.path.isfile(minion_run.summary_file[0]) and not os.path.isdir(minion_run.nanoseq_dir):
        logger.info('Sequencing done for run {}. Attempting to start analysis.'.format(minion_run.run_dir))
        if not minion_run.nanoseq_sample_sheet:
            minion_run.parse_lims_sample_sheet()

        if os.path.isfile(minion_run.nanoseq_sample_sheet):
            if nanoseq_ongoing:
                logger.warn('Nanoseq already started, will not attempt to start for {}'.format(minion_run.run_dir))
            elif sequencing_ongoing:
                logger.warn('Sequencing ongoing, will not attempt to start Nanoseq for {}'.format(minion_run.run_dir))
            else:
                minion_run.start_nanoseq()
                nanoseq_ongoing = True

        else:
            logger.warn('Samplesheet not found for run {}. Operator notified. Skipping.'.format(minion_run.run_dir))
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
                    minion_run.anglerfish_sample_sheet = os.path.join(minion_run.run_dir, 'anglerfish_sample_sheet.csv')  # For cronjob, AF sample sheet was generated at previous run. Consider parsing AF sample sheet separately here instead?
                if os.path.isfile(minion_run.anglerfish_sample_sheet):
                    minion_run.start_anglerfish()
                else:
                    logger.warn('Anglerfish sample sheet missing for run {}. '
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
                            logger.warn('An error occurred during transfer of run {} '
                                        'to Irma. Notifying operator.'.format(minion_run.run_id))
                            email_subject = ('Run processed with errors: {}'.format(minion_run.run_id))
                            email_message = ('Run {} has been analysed, but an error occurred during '
                                             'transfer.').format(minion_run.run_id)
                            send_mail(email_subject, email_message, email_recipients)

                    else:
                        logger.warn('The following run has already been transferred, '
                                    'skipping: {}'.format(minion_run.run_id))

                else:
                    logger.warn('Anglerfish exited with a non-zero exit status for run {}. '
                                'Notifying operator.'.format(minion_run.run_id))
                    email_subject = ('Run processed with errors: {}'.format(minion_run.run_id))
                    email_message = ('Anglerfish exited with errors for run {}. Please '
                                     'check the log files and restart.').format(minion_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)

        else:
            logger.warn('Nanoseq exited with a non-zero exit status for run {}. '
                        'Notifying operator.'.format(minion_run.run_id))
            email_subject = ('Analysis failed for run {}'.format(minion_run.run_id))
            email_message = 'The nanoseq analysis failed for run {}.'.format(minion_run.run_id)
            send_mail(email_subject, email_message, email_recipients)

    else:
        logger.info('Run {} not finished sequencing yet. Skipping.'.format(minion_run.run_id))

    return nanoseq_ongoing

def process_promethion_run(promethion_run):
    """Process promethion runs."""
    email_recipients = CONFIG.get('mail').get('recipients')
    logger.info('Processing run {}'.format(promethion_run.run_id))
    
    if len(promethion_run.summary_file) and os.path.isfile(promethion_run.summary_file[0]):
        logger.info('Sequencing done for run {}. Attempting to start processing.'.format(promethion_run.run_id))
        if promethion_run.is_not_transferred():
            if promethion_run.transfer_run():
                if promethion_run.update_transfer_log():
                    logger.info('Run {} has been synced to the analysis cluster.'.format(promethion_run.run_id))
                else:
                    email_subject = ('Run processed with errors: {}'.format(promethion_run.run_id))
                    email_message = ('Run {} has been transferred, but an error occurred while updating '
                                        'the transfer log').format(promethion_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)
            
                if promethion_run.archive_run():
                    logger.info('Run {} is finished and has been archived. '
                                'Notifying operator.'.format(promethion_run.run_id))
                    email_subject = ('Run successfully processed: {}'.format(promethion_run.run_id))
                    email_message = ('Run {} has been transferred and archived '
                                    'successfully.').format(promethion_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)
                else:
                    email_subject = ('Run processed with errors: {}'.format(promethion_run.run_id))
                    email_message = ('Run {} has been analysed, but an error occurred during '
                                    'archiving').format(promethion_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)

            else:
                email_subject = ('Run processed with errors: {}'.format(promethion_run.run_id))
                email_message = ('An error occurred during transfer of run {} '
                                'to the analysis cluster.').format(promethion_run.run_id)
                send_mail(email_subject, email_message, email_recipients)

        else:
            logger.warn('The following run has already been transferred, '
                        'skipping: {}'.format(promethion_run.run_id))
    else:
        logger.info('Run {} not finished sequencing yet. Skipping.'.format(promethion_run.run_id))


def process_minion_runs(run, nanoseq_sample_sheet, anglerfish_sample_sheet):
    """Find MinION QC runs and kick off processing."""
    if run:
        minion_run = MinION(os.path.abspath(run), nanoseq_sample_sheet, anglerfish_sample_sheet)
        process_minion_run(minion_run)
    else:
        runs_to_process = find_runs_to_process()
        sequencing_ongoing, nanoseq_ongoing = False, False
        for run_dir in runs_to_process:
            if check_ongoing_sequencing(run_dir):
                sequencing_ongoing = True
                break
        for run_dir in runs_to_process:
            minion_run = MinION(run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet)
            nanoseq_ongoing = process_minion_run(minion_run,
                                                 sequencing_ongoing=sequencing_ongoing,
                                                 nanoseq_ongoing=nanoseq_ongoing)

def process_promethion_runs(run):
    """Find promethion runs and kick off processing."""
    if run:
        promethion_run = PromethION(os.path.abspath(run))
        process_promethion_run(promethion_run)
    else:
        # Locate all runs in /srv/ngi_data/sequencing/promethion
        runs_to_process = find_runs_to_process()
        for run_dir in runs_to_process:
            promethion_run = PromethION(run_dir)
            process_promethion_run(promethion_run)