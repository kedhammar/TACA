"""Nanopore analysis methods for TACA."""
import os
import logging

from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
from taca.nanopore import ONT_run

logger = logging.getLogger(__name__)


def find_ont_runs(dir_to_search, skip_dirs):
    """Takes an input dir, expected to contain ONT run dirs.
    Append all found ONT run dirs to a list and return it"""

    try:
        found_run_dirs = []
        for found_dir in os.listdir(dir_to_search):
            if (
                os.path.isdir(os.path.join(dir_to_search, found_dir))
                and found_dir not in skip_dirs
            ):
                found_run_dirs.append(os.path.join(dir_to_search, found_dir))

    except OSError:
        logger.warning(
            "There was an issue locating the following directory: {}. "
            "Please check that it exists and try again.".format(dir_to_search)
        )

    return found_run_dirs


def transfer_ont_run(ont_run):
    """Transfer ONT run to HPC cluster."""

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


def ont_transfer(run_abspath):
    """CLI entry function.

    Find finished ONT runs in ngi-nas and transfer to HPC cluster.
    """

    if run_abspath:
        ont_run = ONT_run(os.path.abspath(run_abspath))
        transfer_ont_run(ont_run)

    else:
        # If no run is specified, locate all runs
        ont_data_dirs = (
            CONFIG.get("nanopore_analysis").get("ont_transfer").get("data_dirs")
        )
        skip_dirs = (
            CONFIG.get("nanopore_analysis").get("ont_transfer").get("ignore_dirs")
        )
        for data_dir in ont_data_dirs:
            runs_to_process = find_ont_runs(data_dir, skip_dirs)
            for run_dir in runs_to_process:
                ont_run = ONT_run(run_dir)
                transfer_ont_run(ont_run)


def ont_updatedb(run_abspath):
    """CLI entry function."""

    ont_run = ONT_run(os.path.abspath(run_abspath))
    ont_run.update_db(force_update=True)