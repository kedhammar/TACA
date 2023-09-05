"""Nanopore analysis methods for TACA."""
import os
import logging
import re
import traceback

from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
from taca.nanopore.ONT_run import ONT_run, ONT_RUN_PATTERN

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
            logger.info(f"Found ONT run {found_dir}.")
            found_run_dirs.append(os.path.join(dir_to_search, found_dir))

    return found_run_dirs


def send_error_mail(run_name, error: BaseException):

    email_subject = f"Run processed with errors: {run_name}"
    email_message = "{}\n\n{}".format(
        str(error),
        traceback.format_exc(),
    )
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

    logger.info(f"{ont_run.run_name}: Touching StatusDB...")
    ont_run.touch_db_entry()
    logger.info(f"{ont_run.run_name}: Touching StatusDB successful...")

    if ont_run.is_synced():
        logger.info(f"{ont_run.run_name}: Run is fully synced.")

        if not ont_run.is_transferred():
            logger.info(f"{ont_run.run_name}: Processing transfer...")

            # Assert all files are in place
            logger.info(f"{ont_run.run_name}: Asserting run contents...")
            ont_run.assert_contents()
            logger.info(f"{ont_run.run_name}: Asserting run contents successful.")

            # Update StatusDB
            logger.info(f"{ont_run.run_name}: Updating StatusDB...")
            ont_run.update_db_entry()
            logger.info(f"{ont_run.run_name}: Updating StatusDB successful.")

            # Transfer HTML report
            logger.info(f"{ont_run.run_name}: Put HTML report on GenStat...")
            ont_run.transfer_html_report()
            logger.info(f"{ont_run.run_name}: Put HTML report on GenStat successful.")

            # Copy metadata
            logger.info(f"{ont_run.run_name}: Copying metadata...")
            ont_run.transfer_metadata()
            logger.info(f"{ont_run.run_name}: Copying metadata successful.")

            # Transfer run
            logger.info(f"{ont_run.run_name}: Transferring to cluster...")
            ont_run.transfer_run()
            logger.info(f"{ont_run.run_name}: Transferring to cluster successful.")

            # Update transfer log
            logger.info(f"{ont_run.run_name}: Updating transfer log...")
            ont_run.update_transfer_log()
            logger.info(f"{ont_run.run_name}: Updating transfer log successful.")

            # Archive run
            logger.info(f"{ont_run.run_name}: Archiving run...")
            ont_run.archive_run()
            logger.info(f"{ont_run.run_name}: Archiving run successful.")

        else:
            logger.warning(
                f"{ont_run.run_name}: Run is already logged as transferred, skipping."
            )
    else:
        logger.info(f"{ont_run.run_name}: Run is not fully synced, skipping.")


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
                    send_error_mail(os.path.basename(run_dir), e)


def ont_updatedb(run_abspath: str):
    """CLI entry function."""

    ont_run = ONT_run(os.path.abspath(run_abspath))

    logger.info(
        f"{ont_run.run_name}: Manually updating StatusDB, ignoring run status..."
    )
    ont_run.update_db(force_update=True)
    logger.info(f"{ont_run.run_name}: Manually updating StatusDB successful.")
