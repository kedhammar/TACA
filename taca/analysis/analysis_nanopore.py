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
    skip_dirs = ['core-dump-db', 'intermediate', 'nosync', 'queued_reads', 'reads', 'user_scripts']
    try:
        found_top_dirs = [os.path.join(nanopore_data_dir, top_dir) for top_dir in os.listdir(nanopore_data_dir)
                            if os.path.isdir(os.path.join(nanopore_data_dir, top_dir))
                            and top_dir not in skip_dirs]
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

def check_ongoing_sequencing(run_dir):
    """Check if sequencing is ongoing for the given run dir """
    summary_file = glob.glob(run_dir + '/final_summary*.txt')

    if not len(summary_file):
        logger.info('Sequencing ongoing for run {}. Will not start nanoseq at this time'.format(run_dir))
        return True
    else:
        return False

def process_minion_run(MinionRun, sequencing_ongoing=False, nanoseq_ongoing=False): 
    """Process minion runs.
    
    Will not start nanoseq if a sequencing run is ongoing, to limit memory usage.
    Will also maximum start one nanoseq run at once, for the same reason.
    """
    qc_run = True
    if MinionRun.nanoseq_sample_sheet and not MinionRun.anglerfish_sample_sheet:
        qc_run = False

    logger.info('Processing run: {} as a {}'.format(MinionRun.run_dir, 'QC run' if qc_run else 'non-QC run'))
    email_recipients = CONFIG.get('mail').get('recipients')

    if len(MinionRun.summary_file) and os.path.isfile(MinionRun.summary_file[0]) and not os.path.isdir(MinionRun.nanoseq_dir):
        logger.info('Sequencing done for run {}. Attempting to start analysis.'.format(MinionRun.run_dir))
        if not MinionRun.nanoseq_sample_sheet:
            MinionRun.parse_lims_sample_sheet()

        if os.path.isfile(MinionRun.nanoseq_sample_sheet):
            if nanoseq_ongoing:
                logger.warn('Nanoseq already started, will not attempt to start for {}'.format(MinionRun.run_dir))
            elif sequencing_ongoing:
                logger.warn('Sequencing ongoing, will not attempt to start Nanoseq for {}'.format(MinionRun.run_dir))
            else:
                MinionRun.start_nanoseq()
                nanoseq_ongoing = True

        else:
            logger.warn('Samplesheet not found for run {}. Operator notified. Skipping.'.format(MinionRun.run_dir))
            email_subject = ('Samplesheet missing for run {}'.format(os.path.basename(MinionRun.run_dir)))
            email_message = 'There was an issue locating the samplesheet for run {}.'.format(MinionRun.run_dir)
            send_mail(email_subject, email_message, email_recipients)

    elif os.path.isdir(MinionRun.nanoseq_dir) and not os.path.isfile(MinionRun.nanoseq_exit_status_file):
        logger.info('Nanoseq has started for run {} but is not yet done. Skipping.'.format(MinionRun.run_dir))

    elif os.path.isdir(MinionRun.nanoseq_dir) and os.path.isfile(MinionRun.nanoseq_exit_status_file):
        nanoseq_successful = MinionRun.check_exit_status(MinionRun.nanoseq_exit_status_file)
        if nanoseq_successful:
            if qc_run and not os.path.isdir(MinionRun.anglerfish_dir):
                logger.info('Nanoseq done for run {}. Attempting to start Anglerfish.'.format(MinionRun.run_id))
                if not MinionRun.anglerfish_sample_sheet:
                    MinionRun.anglerfish_sample_sheet = os.path.join(MinionRun.run_dir, 'anglerfish_sample_sheet.csv')
                MinionRun.start_anglerfish()

            elif qc_run and not os.path.isfile(MinionRun.anglerfish_exit_status_file):
                logger.info('Anglerfish has started for run {} but is not yet done. Skipping.'.format(MinionRun.run_id))

            elif qc_run and os.path.isfile(MinionRun.anglerfish_exit_status_file):
                anglerfish_successful = MinionRun.check_exit_status(MinionRun.anglerfish_exit_status_file)
                if anglerfish_successful:
                    MinionRun.copy_results_for_lims()
                    logger.info('Anglerfish finished OK for run {}. Notifying operator.'.format(MinionRun.run_id))
                    email_subject = ('Anglerfish successfully processed run {}'.format(MinionRun.run_id))
                    email_message = ('Anglerfish has successfully finished for run {}. Please '
                                     'finish the QC step in lims.').format(MinionRun.run_id)
                    send_mail(email_subject, email_message, email_recipients)

                    if MinionRun.is_not_transferred():
                        if MinionRun.transfer_run():
                            MinionRun.update_transfer_log()
                            logger.info('Run {} has been synced to the analysis cluster.'.format(MinionRun.run_id))
                            MinionRun.archive_run()
                            logger.info('Run {} is finished and has been archived. Notifying operator.'.format(MinionRun.run_id))
                            email_subject = ('Run successfully processed: {}'.format(MinionRun.run_id))
                            email_message = ('Run {} has been analysed, transferred and archived '
                                             'successfully.').format(MinionRun.run_id)
                            send_mail(email_subject, email_message, email_recipients)

                        else:
                            logger.warn('An error occurred during transfer of run {} '
                                        'to Irma. Notifying operator.'.format(MinionRun.run_id))
                            email_subject = ('Run processed with errors: {}'.format(MinionRun.run_id))
                            email_message = ('Run {} has been analysed, but an error occurred during '
                                             'transfer.').format(MinionRun.run_id)
                            send_mail(email_subject, email_message, email_recipients)

                    else:
                        logger.warn('The following run has already been transferred, '
                                    'skipping: {}'.format(MinionRun.run_id))

                else:
                    logger.warn('Anglerfish exited with a non-zero exit status for run {}. '
                                'Notifying operator.'.format(MinionRun.run_id))
                    email_subject = ('Run processed with errors: {}'.format(MinionRun.run_id))
                    email_message = ('Anglerfish exited with errors for run {}. Please '
                                     'check the log files and restart.').format(MinionRun.run_id)
                    send_mail(email_subject, email_message, email_recipients)

            elif not qc_run:
                if MinionRun.is_not_transferred():
                    if MinionRun.transfer_run():
                        MinionRun.update_transfer_log()
                        logger.info('Run {} has been synced to the analysis cluster.'.format(MinionRun.run_id))
                        MinionRun.archive_run()
                        logger.info('Run {} is finished and has been archived. Notifying operator.'.format(MinionRun.run_id))
                        email_subject = ('Run successfully processed: {}'.format(MinionRun.run_id))
                        email_message = ('Run {} has been analysed, transferred and archived '
                                         'successfully.').format(MinionRun.run_id)
                        send_mail(email_subject, email_message, email_recipients)

                    else:
                        logger.warn('An error occurred during transfer of run {} '
                                    'to Irma. Notifying operator.'.format(MinionRun.run_id))
                        email_subject = ('Run processed with errors: {}'.format(MinionRun.run_id))
                        email_message = ('Run {} has been analysed, but an error occurred during '
                                         'transfer.').format(MinionRun.run_id)
                        send_mail(email_subject, email_message, email_recipients)

                else:
                    logger.warn('The following run has already been transferred, '
                                'skipping: {}'.format(MinionRun.run_id))

        else:
            logger.warn('Nanoseq exited with a non-zero exit status for run {}. '
                        'Notifying operator.'.format(MinionRun.run_id))
            email_subject = ('Analysis failed for run {}'.format(MinionRun.run_id))
            email_message = 'The nanoseq analysis failed for run {}.'.format(MinionRun.run_id)
            send_mail(email_subject, email_message, email_recipients)

    else:
        logger.info('Run {} not finished sequencing yet. Skipping.'.format(MinionRun.run_id))

    return nanoseq_ongoing

def process_promethion_run(promethion_run):
    """Process promethion runs."""
    # WIP - please ignore
    email_recipients = CONFIG.get('mail').get('recipients')
    logger.info('Processing run {}'.format(promethion_run.run_id))
    
    if len(promethion_run.summary_file) and os.path.isfile(promethion_run.summary_file[0]):
        logger.info('Sequencing done for run {}. Attempting to start processing.'.format(promethion_run.run_id))
        if promethion_run.is_not_transferred():
            if promethion_run.transfer_run():
                promethion_run.update_transfer_log()
                logger.info('Run {} has been synced to the analysis cluster.'.format(promethion_run.run_id))
            
                promethion_run.archive_run()
                logger.info('Run {} is finished and has been archived. '
                            'Notifying operator.'.format(promethion_run.run_id))
                email_subject = ('Run successfully processed: {}'.format(promethion_run.run_id))
                email_message = ('Run {} has been transferred and archived '
                                'successfully.').format(promethion_run.run_id)
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
    """Find minion runs and kick off processing."""
    if run:
        MinionRun = MinION(os.path.abspath(run), nanoseq_sample_sheet, anglerfish_sample_sheet)
        process_minion_run(MinionRun)
    else:
        runs_to_process = find_runs_to_process()
        sequencing_ongoing, nanoseq_ongoing = False, False
        for run_dir in runs_to_process:
            if check_ongoing_sequencing(run_dir):
                sequencing_ongoing = True
                break
        for run_dir in runs_to_process:
            MinionRun = MinION(run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet)
            nanoseq_ongoing = process_minion_run(MinionRun,
                                                 sequencing_ongoing=sequencing_ongoing, 
                                                 nanoseq_ongoing=nanoseq_ongoing)

def process_promethion_runs(run):
    """Find promethion runs and kick off processing."""
    if run:
        promethion_run = PromethION(os.path.abspath(run))
        process_promethion_run(promethion_run)
    else:
        # Locate all runs in /srv/ngi_data
        runs_to_process = find_runs_to_process()
        for run_dir in runs_to_process:
            promethion_run = PromethION(run_dir)
            process_promethion_run(promethion_run)