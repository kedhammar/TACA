"""Nanopore analysis methods for TACA."""
import os
import logging
import re
import traceback

from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
from taca.nanopore import ONT_run, ONT_RUN_PATTERN

logger = logging.getLogger(__name__)


def find_run_dirs(dir_to_search: str, skip_dirs: list):
    """Takes an input dir, expected to contain ONT run dirs.
    Append all found ONT run dirs to a list and return it"""

    found_run_dirs = []
    for found_dir in os.listdir(dir_to_search):
        if (
            os.path.isdir(os.path.join(dir_to_search, found_dir))
            and found_dir not in skip_dirs
            and re.match(ONT_RUN_PATTERN, found_dir)
        ):
            found_run_dirs.append(os.path.join(dir_to_search, found_dir))

    return found_run_dirs


def send_error_mail(ont_run: ONT_run, error: BaseException):

    email_subject = f"Run processed with errors: {ont_run.run_id}"
    email_message = f"{str(error)}\n\n{traceback.format_exc(error)}"
    email_recipients = CONFIG.get("mail").get("recipients")

    send_mail(email_subject, email_message, email_recipients)


def transfer_ont_run(ont_run: ONT_run):
    """Transfer ONT run to HPC cluster."""

    logger.info("Processing run {}".format(ont_run.run_id))

    # Update StatusDB
    ont_run.update_db()

    if ont_run.is_synced:
        logger.info(f"Run {ont_run.run_id} has finished sequencing.")

        if not ont_run.is_transferred():
            logger.info(f"Run {ont_run.run_id} is not yet transferred.")

            # Copy metadata
            logger.info(f"Copying metadata of {ont_run.run_id}.")
            ont_run.transfer_metadata()
            logger.info(
                f"Metadata of run {ont_run.run_id} has been synced to {ont_run.metadata_dir}"
            )

            # Transfer run
            if ont_run.transfer_run():
                if ont_run.update_transfer_log():
                    logger.info(
                        f"Run {ont_run.run_id} has been synced to the analysis cluster."
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


def ont_transfer(run_abspath: str or None):
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
            run_dirs = find_run_dirs(data_dir, skip_dirs)
            for run_dir in run_dirs:
                # Send error mails at run-level
                try:
                    ont_run = ONT_run(run_dir)
                    transfer_ont_run(ont_run)
                except BaseException as e:
                    send_error_mail(ont_run, e)


def ont_updatedb(run_abspath: str):
    """CLI entry function."""

    ont_run = ONT_run(os.path.abspath(run_abspath))
    ont_run.update_db(force_update=True)