"""Nanopore analysis methods for TACA."""
import os
import logging
import re
import traceback

from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
from taca.nanopore.ONT_run_classes import (
    ONT_run,
    ONT_user_run,
    ONT_qc_run,
    ONT_RUN_PATTERN,
)

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
            logger.info(f"Found ONT run {found_dir} in {dir_to_search}")
            found_run_dirs.append(os.path.join(dir_to_search, found_dir))

    return found_run_dirs


def send_error_mail(run_name, error: BaseException):

    email_subject = f"Run processed with errors: {run_name}"
    email_message = "{}\n\n{}".format(
        str(error),
        traceback.format_exc(),
    )
    email_recipients = CONFIG["mail"]["recipients"]

    send_mail(email_subject, email_message, email_recipients)


def process_user_run(ONT_user_run: ONT_run):
    """This control function orchestrates the sequential execution of the ONT_user_run class methods.

    For a single ONT user run...

        a) If not finished:
            - Ensure there is a database entry corresponding to an ongoing run

        b) If finished:
            - Ensure all necessary files to proceed with processing are present
            - Ensure there is a database entry corresponding to an ongoing run
            - Update the StatusDB entry
            - Copy metadata
            - Copy HTML report to GenStat
            - Transfer run to cluster
            - Update transfer log
            - Archive run

    Any errors raised here-in should be sent with traceback as an email.
    """

    logger.info(f"{ONT_user_run.run_name}: Touching StatusDB...")
    ONT_user_run.touch_db_entry()
    logger.info(f"{ONT_user_run.run_name}: Touching StatusDB successful...")

    if ONT_user_run.is_synced():
        logger.info(f"{ONT_user_run.run_name}: Run is fully synced.")

        if not ONT_user_run.is_transferred():
            logger.info(f"{ONT_user_run.run_name}: Processing transfer...")

            # Assert all files are in place
            logger.info(f"{ONT_user_run.run_name}: Asserting run contents...")
            ONT_user_run.assert_contents()
            logger.info(f"{ONT_user_run.run_name}: Asserting run contents successful.")

            # Update StatusDB
            logger.info(f"{ONT_user_run.run_name}: Updating StatusDB...")
            ONT_user_run.update_db_entry()
            logger.info(f"{ONT_user_run.run_name}: Updating StatusDB successful.")

            # Copy HTML report
            logger.info(f"{ONT_user_run.run_name}: Put HTML report on GenStat...")
            ONT_user_run.copy_html_report()
            logger.info(
                f"{ONT_user_run.run_name}: Put HTML report on GenStat successful."
            )

            # Copy metadata
            logger.info(f"{ONT_user_run.run_name}: Copying metadata...")
            ONT_user_run.copy_metadata()
            logger.info(f"{ONT_user_run.run_name}: Copying metadata successful.")

            # Transfer run
            logger.info(f"{ONT_user_run.run_name}: Transferring to cluster...")
            ONT_user_run.transfer_run()
            logger.info(f"{ONT_user_run.run_name}: Transferring to cluster successful.")

            # Update transfer log
            logger.info(f"{ONT_user_run.run_name}: Updating transfer log...")
            ONT_user_run.update_transfer_log()
            logger.info(f"{ONT_user_run.run_name}: Updating transfer log successful.")

            # Archive run
            logger.info(f"{ONT_user_run.run_name}: Archiving run...")
            ONT_user_run.archive_run()
            logger.info(f"{ONT_user_run.run_name}: Archiving run successful.")

        else:
            logger.warning(
                f"{ONT_user_run.run_name}: Run is already logged as transferred, skipping."
            )
    else:
        logger.info(f"{ONT_user_run.run_name}: Run is not fully synced, skipping.")


def process_qc_run(ont_qc_run: ONT_qc_run):
    """This control function orchestrates the sequential execution of the ONT_qc_run class methods.

    Any errors raised here-in should be sent with traceback as an email.
    """

    logger.info(f"{ont_qc_run.run_name}: Touching StatusDB...")
    ont_qc_run.touch_db_entry()
    logger.info(f"{ont_qc_run.run_name}: Touching StatusDB successful...")

    if ont_qc_run.is_synced():
        logger.info(f"{ont_qc_run.run_name}: Run is fully synced.")

        logger.info(
            f"{ont_qc_run.run_name}: Run is QC, checking for Anglerfish samplesheet..."
        )

        if ont_qc_run.fetch_anglerfish_samplesheet():
            logger.info(
                f"{ont_qc_run.run_name}: Anglerfish samplesheet found and copied to run dir."
            )

            logger.info(f"{ont_qc_run.run_name}: Running Anglerfish...")
            ont_qc_run.run_anglerfish()
            logger.info(f"{ont_qc_run.run_name}: Running Anglerfish successful.")

        else:
            logger.info(
                f"{ont_qc_run.run_name}: Anglerfish samplesheet not found, skipping."
            )


def process_run(run_abspath: str):
    """This gate function instantiates an appropriate subclass
    for the ONT run and processes it"""

    ont_run = ONT_run(run_abspath)

    if ont_run.qc:
        process_qc_run(ONT_qc_run(run_abspath))
    else:
        process_user_run(ONT_user_run(run_abspath))


def ont_transfer(run_abspath: str or None):
    """CLI entry function.

    Find finished ONT runs in ngi-nas and transfer to HPC cluster.
    """

    if run_abspath:
        process_run(run_abspath)

    # If no run is specified, locate all runs
    else:

        for run_type in ["user_run", "qc_run"]:

            logger.info(f"Looking for runs of type '{run_type}'...")

            data_dirs = CONFIG["nanopore_analysis"][run_type]["data_dirs"]
            ignore_dirs = CONFIG["nanopore_analysis"][run_type]["ignore_dirs"]

            for data_dir in data_dirs:

                run_dirs = find_run_dirs(data_dir, ignore_dirs)
                for run_dir in run_dirs:

                    # Send error mails at run-level
                    try:
                        process_run(run_dir)

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
