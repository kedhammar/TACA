"""Nanopore analysis methods for TACA."""
import os
import logging
import subprocess

from dateutil.parser import parse
from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
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


def find_ont_transfer_runs(ont_data_dir, skip_dirs):
    """Find runs in ngi-nas.
    These are assumed to be flowcell dirs, not project dirs.
    """
    try:
        found_dirs = [
            os.path.join(ont_data_dir, top_dir)
            for top_dir in os.listdir(ont_data_dir)
            if os.path.isdir(os.path.join(ont_data_dir, top_dir))
            and top_dir not in skip_dirs
        ]
    except OSError:
        logger.warning(
            "There was an issue locating the following directory: {}. "
            "Please check that it exists and try again.".format(ont_data_dir)
        )
    return found_dirs


def process_minion_qc_run(minion_run):
    """Process MinION QC runs on Squiggle.
    """
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
            
        if minion_run.anglerfish_sample_sheet and os.path.isfile(minion_run.anglerfish_sample_sheet):
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
                email_subject = (
                    "Anglerfish successfully processed run {}".format(
                        minion_run.run_id
                    )
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
                        send_mail(
                            email_subject, email_message, email_recipients
                        )

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
                        send_mail(
                            email_subject, email_message, email_recipients
                        )
                    else:
                        email_subject = "Run processed with errors: {}".format(
                            minion_run.run_id
                        )
                        email_message = (
                            "Run {} has been analysed, but an error occurred during "
                            "archiving"
                        ).format(minion_run.run_id)
                        send_mail(
                            email_subject, email_message, email_recipients
                        )
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
            email_subject = "Run processed with errors: {}".format(
                minion_run.run_id
            )
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


def transfer_ont_run(ont_run):
    """Transfer ONT runs to HPC cluster."""
    email_recipients = CONFIG.get("mail").get("recipients")
    logger.info("Processing run {}".format(ont_run.run_id))

    # Update StatusDB
    try:
        ont_run.update_db()
    except Exception as e:
        logger.warning(f"Database update for run {ont_run.run_id} failed")
        email_subject = "Run processed with errors: {}".format(ont_run.run_id)
        email_message = (
            f"An error occured when updating statusdb with run {ont_run.run_id}.\n{e}"
        )
        send_mail(email_subject, email_message, email_recipients)

    if os.path.isfile(ont_run.sync_finished_indicator):
        logger.info(
            "Sequencing done for run {}. Attempting to start processing.".format(
                ont_run.run_id
            )
        )
        if ont_run.is_not_transferred():

            # Copy metadata
            try:
                ont_run.transfer_metadata()
                logger.info(
                    f"Metadata of run {ont_run.run_id} has been synced to {ont_run.metadata_dir}"
                )
            except BaseException as e:
                email_subject = f"Run processed with errors: {ont_run.run_id}"
                email_message = f"Run {ont_run.run_id} has been analysed, but an error occurred when copying the metadata: \n{str(e)}"
                send_mail(email_subject, email_message, email_recipients)

            # Transfer run
            if ont_run.transfer_run():
                if ont_run.update_transfer_log():
                    logger.info(
                        "Run {} has been synced to the analysis cluster.".format(
                            ont_run.run_id
                        )
                    )
                else:
                    email_subject = "Run processed with errors: {}".format(
                        ont_run.run_id
                    )
                    email_message = (
                        "Run {} has been transferred, but an error occurred while updating "
                        "the transfer log"
                    ).format(ont_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)

                if ont_run.archive_run():
                    logger.info(
                        "Run {} is finished and has been archived. "
                        "Notifying operator.".format(ont_run.run_id)
                    )
                    email_subject = "Run successfully processed: {}".format(
                        ont_run.run_id
                    )
                    email_message = (
                        "Run {} has been transferred and archived " "successfully."
                    ).format(ont_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)
                else:
                    email_subject = "Run processed with errors: {}".format(
                        ont_run.run_id
                    )
                    email_message = (
                        "Run {} has been analysed, but an error occurred during "
                        "archiving"
                    ).format(ont_run.run_id)
                    send_mail(email_subject, email_message, email_recipients)

            else:
                email_subject = "Run processed with errors: {}".format(ont_run.run_id)
                email_message = (
                    "An error occurred during transfer of run {} "
                    "to the analysis cluster."
                ).format(ont_run.run_id)
                send_mail(email_subject, email_message, email_recipients)

        else:
            logger.warning(
                "The following run has already been transferred, "
                "skipping: {}".format(ont_run.run_id)
            )
    else:
        logger.info(
            "Run {} not finished sequencing yet. Skipping.".format(ont_run.run_id)
        )


def process_minion_qc_runs(run, anglerfish_sample_sheet):
    """Find and process MinION QC runs on Squiggle."""
    if run:
        if is_date(os.path.basename(run).split("_")[0]):
            minion_run = MinIONqc(
                os.path.abspath(run), anglerfish_sample_sheet
            )
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
            minion_run = MinIONqc(
                run_dir, anglerfish_sample_sheet
            )
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


def transfer_finished(run):
    """Find finished ONT runs in ngi-nas and transfer to HPC cluster."""
    if run:
        if is_date(os.path.basename(run).split("_")[0]):
            if "minion" in run:
                ont_run = MinionTransfer(os.path.abspath(run))
            elif "promethion" in run:
                ont_run = PromethionTransfer(os.path.abspath(run))
            transfer_ont_run(ont_run)
        else:
            logger.warning(
                "The specified path is not a flow cell. Please "
                "provide the full path to the flow cell you wish to process."
            )
    else:
        # Locate all runs in /srv/ngi_data/sequencing/promethion and /srv/ngi_data/sequencing/minion
        ont_data_dirs = (
            CONFIG.get("nanopore_analysis").get("ont_transfer").get("data_dirs")
        )
        skip_dirs = (
            CONFIG.get("nanopore_analysis").get("ont_transfer").get("ignore_dirs")
        )
        for data_dir in ont_data_dirs:
            runs_to_process = find_ont_transfer_runs(data_dir, skip_dirs)
            for run_dir in runs_to_process:
                if "minion" in data_dir:
                    ont_run = MinionTransfer(run_dir)
                    transfer_ont_run(ont_run)
                elif "promethion" in data_dir:
                    ont_run = PromethionTransfer(run_dir)
                    transfer_ont_run(ont_run)


def ont_updatedb_from_cli(run):

    if is_date(os.path.basename(run).split("_")[0]):
        if "minion" in run:
            ont_run = MinionTransfer(os.path.abspath(run))
        elif "promethion" in run:
            ont_run = PromethionTransfer(os.path.abspath(run))
        ont_run.update_db(force_update=True)
    else:
        logger.warning(
            "The specified path is not a flow cell. Please "
            "provide the full path to the flow cell you wish to process."
        )
