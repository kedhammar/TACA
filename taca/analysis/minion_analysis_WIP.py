import subprocess
import os
import logging

from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
from taca.nanopore.minion_run_class import MinIONdelivery, MinIONqc


logger = logging.getLogger(__name__)


def find_minion_runs(minion_data_dir, skip_dirs):
    """Find nanopore runs to process."""
    found_run_dirs = []
    try:
        found_top_dirs = [
            os.path.join(minion_data_dir, top_dir)
            for top_dir in os.listdir(minion_data_dir)
            if os.path.isdir(os.path.join(minion_data_dir, top_dir))
            and top_dir not in skip_dirs
        ]
    except OSError:
        logger.warning(
            "There was an issue locating the following directory: {}. "
            "Please check that it exists and try again.".format(minion_data_dir)
        )
    # Get the actual location of the run directories in /var/lib/MinKnow/data/QC_runs/USERDETERMINEDNAME/USERDETSAMPLENAME/run
    if found_top_dirs:
        for top_dir in found_top_dirs:
            if os.path.isdir(top_dir):
                for sample_dir in os.listdir(top_dir):
                    if os.path.isdir(os.path.join(top_dir, sample_dir)):
                        for run_dir in os.listdir(os.path.join(top_dir, sample_dir)):
                            found_run_dirs.append(
                                os.path.join(top_dir, sample_dir, run_dir)
                            )
    else:
        logger.warning(
            "Could not find any run directories in {}".format(minion_data_dir)
        )
    return found_run_dirs


def process_minion_qc_run(minion_run):
    """Process MinION QC runs on Squiggle."""
    logger.info("Processing QC run: {}".format(minion_run.run_dir))
    email_recipients = CONFIG.get("mail").get("recipients")
    if not len(minion_run.summary_file):
        logger.info(
            "Sequencing is still ongoing for run {}. Skipping.".format(
                minion_run.run_id
            )
        )
        return

    if (
        len(minion_run.summary_file)
        and os.path.isfile(minion_run.summary_file[0])
        and not os.path.isdir(minion_run.anglerfish_dir)
    ):
        logger.info(
            "Sequencing is done for run {}. Attempting to start Anglerfish.".format(
                minion_run.run_id
            )
        )
        if not minion_run.anglerfish_sample_sheet:
            minion_run.anglerfish_sample_sheet = minion_run.get_anglerfish_samplesheet()

        if minion_run.anglerfish_sample_sheet and os.path.isfile(
            minion_run.anglerfish_sample_sheet
        ):
            minion_run.start_anglerfish()
        else:
            logger.warning(
                "Anglerfish sample sheet missing for run {}. "
                "Please provide one using --anglerfish_sample_sheet "
                "or complete the correct lims step.".format(minion_run.run_id)
            )
    elif not os.path.isfile(minion_run.anglerfish_exit_status_file):
        logger.info(
            "Anglerfish has started for run {} but is not yet done. Skipping.".format(
                minion_run.run_id
            )
        )
    elif os.path.isfile(minion_run.anglerfish_exit_status_file):
        anglerfish_successful = minion_run.check_exit_status(
            minion_run.anglerfish_exit_status_file
        )
        if anglerfish_successful:
            if minion_run.copy_results_for_lims():
                logger.info(
                    "Anglerfish finished OK for run {}. Notifying operator.".format(
                        minion_run.run_id
                    )
                )
                email_subject = "Anglerfish successfully processed run {}".format(
                    minion_run.run_id
                )
                email_message = (
                    "Anglerfish has successfully finished for run {}. Please "
                    "finish the QC step in lims."
                ).format(minion_run.run_id)
                send_mail(email_subject, email_message, email_recipients)
            else:
                email_subject = "Run processed with errors: {}".format(
                    minion_run.run_id
                )
                email_message = (
                    "Anglerfish has successfully finished for run {} but an error "
                    "occurred while transferring the results to lims."
                ).format(minion_run.run_id)
                send_mail(email_subject, email_message, email_recipients)

            if minion_run.is_not_transferred():
                if minion_run.transfer_run():
                    if minion_run.update_transfer_log():
                        logger.info(
                            "Run {} has been synced to the analysis cluster.".format(
                                minion_run.run_id
                            )
                        )
                    else:
                        email_subject = "Run processed with errors: {}".format(
                            minion_run.run_id
                        )
                        email_message = (
                            "Run {} has been transferred, but an error occurred while updating "
                            "the transfer log"
                        ).format(minion_run.run_id)
                        send_mail(email_subject, email_message, email_recipients)

                    if minion_run.archive_run():
                        logger.info(
                            "Run {} is finished and has been archived. Notifying operator.".format(
                                minion_run.run_id
                            )
                        )
                        email_subject = "Run successfully processed: {}".format(
                            minion_run.run_id
                        )
                        email_message = (
                            "Run {} has been analysed, transferred and archived "
                            "successfully."
                        ).format(minion_run.run_id)
                        send_mail(email_subject, email_message, email_recipients)
                    else:
                        email_subject = "Run processed with errors: {}".format(
                            minion_run.run_id
                        )
                        email_message = (
                            "Run {} has been analysed, but an error occurred during "
                            "archiving"
                        ).format(minion_run.run_id)
                        send_mail(email_subject, email_message, email_recipients)
                else:
                    logger.warning(
                        "An error occurred during transfer of run {} "
                        "to the analysis cluster. Notifying operator.".format(
                            minion_run.run_id
                        )
                    )
                    email_subject = "Run processed with errors: {}".format(
                        minion_run.run_id
                    )
                    email_message = (
                        "Run {} has been analysed, but an error occurred during "
                        "transfer to the analysis cluster."
                    ).format(minion_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)
            else:
                logger.warning(
                    "The following run has already been transferred, "
                    "skipping: {}".format(minion_run.run_id)
                )

        else:
            logger.warning(
                "Anglerfish exited with a non-zero exit status for run {}. "
                "Notifying operator.".format(minion_run.run_id)
            )
            email_subject = "Run processed with errors: {}".format(minion_run.run_id)
            email_message = (
                "Anglerfish exited with errors for run {}. Please "
                "check the log files and restart."
            ).format(minion_run.run_id)
            send_mail(email_subject, email_message, email_recipients)

    return


def process_minion_delivery_run(minion_run):
    """Process minion delivery runs on Squiggle."""
    email_recipients = CONFIG.get("mail").get("recipients")
    logger.info("Processing run {}".format(minion_run.run_id))
    minion_run.dump_path()
    if not len(minion_run.summary_file):  # Run not finished, only rsync
        minion_run.transfer_run()
    else:  # Run finished, rsync and archive
        if minion_run.transfer_run():
            finished_indicator = minion_run.write_finished_indicator()
            destination = os.path.join(
                minion_run.transfer_details.get("destination"), minion_run.run_id
            )
            sync_finished_indicator = ["rsync", finished_indicator, destination]
            process_handle = subprocess.run(sync_finished_indicator)
            minion_run.archive_run()
            logger.info("Run {} has been fully transferred.".format(minion_run.run_id))
            email_subject = "Run successfully processed: {}".format(minion_run.run_id)
            email_message = (
                "Run {} has been transferred and archived " "successfully."
            ).format(minion_run.run_id)
            send_mail(email_subject, email_message, email_recipients)
        else:
            logger.warning(
                "An error occurred during transfer of run {}.".format(minion_run.run_id)
            )
            email_subject = "Run processed with errors: {}".format(minion_run.run_id)
            email_message = (
                "An error occurred during the " "transfer of run {}."
            ).format(minion_run.run_id)
            send_mail(email_subject, email_message, email_recipients)


def process_minion_qc_runs(run, anglerfish_sample_sheet):
    """Find and process MinION QC runs on Squiggle."""
    if run:
        if is_date(os.path.basename(run).split("_")[0]):
            minion_run = MinIONqc(os.path.abspath(run), anglerfish_sample_sheet)
            process_minion_qc_run(minion_run)
        else:
            logger.warning(
                "The specified path is not a flow cell. Please "
                "provide the full path to the flow cell you wish to process."
            )
    else:
        nanopore_data_dir = (
            CONFIG.get("nanopore_analysis").get("minion_qc_run").get("data_dir")
        )
        skip_dirs = (
            CONFIG.get("nanopore_analysis").get("minion_qc_run").get("ignore_dirs")
        )
        runs_to_process = find_minion_runs(nanopore_data_dir, skip_dirs)
        for run_dir in runs_to_process:
            minion_run = MinIONqc(run_dir, anglerfish_sample_sheet)
            process_minion_qc_run(minion_run)


def process_minion_delivery_runs(run):
    """Find MinION delivery runs on Squiggle and transfer them to ngi-nas."""
    if run:
        if is_date(os.path.basename(run).split("_")[0]):
            minion_run = MinIONdelivery(os.path.abspath(run))
            process_minion_delivery_run(minion_run)
        else:
            logger.warning(
                "The specified path is not a flow cell. Please "
                "provide the full path to the flow cell you wish to process."
            )
    else:
        minion_data_dir = (
            CONFIG.get("nanopore_analysis").get("minion_delivery_run").get("data_dir")
        )
        skip_dirs = (
            CONFIG.get("nanopore_analysis")
            .get("minion_delivery_run")
            .get("ignore_dirs")
        )
        runs_to_process = find_minion_runs(minion_data_dir, skip_dirs)
        for run_dir in runs_to_process:
            minion_run = MinIONdelivery(run_dir)
            process_minion_delivery_run(minion_run)
