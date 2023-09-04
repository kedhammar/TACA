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

    logger.info(f"Looking for ONT runs in {dir_to_search}...")

    found_run_dirs = []
    for found_dir in os.listdir(dir_to_search):
        if (
            os.path.isdir(os.path.join(dir_to_search, found_dir))
            and found_dir not in skip_dirs
            and re.match(ONT_RUN_PATTERN, found_dir)
        ):
            logger.info(f"Found ONT run {found_dir} in {dir_to_search}.")
            found_run_dirs.append(os.path.join(dir_to_search, found_dir))

    return found_run_dirs


def send_error_mail(ont_run: ONT_run, error: BaseException):

    email_subject = f"Run processed with errors: {ont_run.run_id}"
    email_message = f"{str(error)}\n\n{traceback.format_exc(error)}"
    email_recipients = CONFIG.get("mail").get("recipients")

    send_mail(email_subject, email_message, email_recipients)


def transfer_ont_run(ont_run: ONT_run):
    """This function orchestrates the sequential execution of the ONT_run class methods.

    For a single ONT run...

        a) If not finished:
            - Ensure there is a database entry corresponding to an ongoing run

        b) If finished:
            - Ensure there is a database entry corresponding to an ongoing run
            - Update the StatusDB entry
            - Copy metadata
            - Copy HTML report to GenStat
            - Transfer run to cluster
            - Update transfer log
            - Archive run

    Any errors raised here-in should be sent with traceback as an email.
    """

    logger.info(f"{ont_run.run_id}: Inspecting StatusDB...")
    ont_run.touch_db_entry()

    if ont_run.is_synced:
        logger.info(f"{ont_run.run_id}: Finished sequencing.")

        if not ont_run.is_transferred():
            logger.info(f"{ont_run.run_id}: Processing...")

            # Update StatusDB
            logger.info(f"{ont_run.run_id}: Updating StatusDB...")
            ont_run.update_db_entry()
            logger.info(f"{ont_run.run_id}: Updating StatusDB successful.")

            # Transfer HTML report
            logger.info(f"{ont_run.run_id}: Copying HTML report to GenStat...")
            ont_run.transfer_html_report()
            logger.info(f"{ont_run.run_id}: Copying HTML report to GenStat successful.")

            # Copy metadata
            logger.info(f"{ont_run.run_id}: Copying metadata...")
            ont_run.transfer_metadata()
            logger.info(f"{ont_run.run_id}: Copying metadata successful.")

            # Transfer run
            logger.info(f"{ont_run.run_id}: Transferring to cluster...")
            ont_run.transfer_run()
            logger.info(f"{ont_run.run_id}: Transferring to cluster successful.")

            # Update transfer log
            logger.info(f"{ont_run.run_id}: Updating transfer log...")
            ont_run.update_transfer_log()
            logger.info(f"{ont_run.run_id}: Updating transfer log successful.")

            # Archive run
            logger.info(f"{ont_run.run_id}: Archiving run...")
            ont_run.update_transfer_log()
            logger.info(f"{ont_run.run_id}: Archiving run successful.")

        else:
            logger.warning(
                f"{ont_run.run_id}: Already logged as transferred, skipping."
            )
    else:
        logger.info(f"{ont_run.run_id}: Not finished sequencing yet, skipping.")


def ont_transfer(run_abspath: str or None):
    """CLI entry function.

    Find finished ONT runs in ngi-nas and transfer to HPC cluster.
    """

    if run_abspath:
        # No need to send mails for manual command, CLI will show errors
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

    logger.info(f"{ont_run.run_id}: Manually updating StatusDB, ignoring run status...")
    ont_run.update_db(force_update=True)
    logger.info(f"{ont_run.run_id}: Manually updating StatusDB successful.")
