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
    email_message = f"{str(error)}\n\n{traceback.format_exc()}"
    email_recipients = CONFIG["mail"]["recipients"]

    logger.warning(
        f"Error encountered for run {run_name}.\nTraceback:\n{email_message}"
    )

    send_mail(email_subject, email_message, email_recipients)

    logger.info("Sent error mail.")


def process_user_run(ont_user_run: ONT_user_run):
    """This control function orchestrates the sequential execution of the ONT_user_run class methods.

    For a single ONT user run...

        - Ensure there is a database entry corresponding to an ongoing run

        If not fully synced:
            - Skip
        If fully synced:
            - Ensure all necessary files to proceed with processing are present
            - Update the StatusDB entry
            - Copy metadata
            - Copy HTML report to GenStat
            - Transfer run to cluster
            - Update transfer log
            - Archive run

    Any errors raised here-in should be sent with traceback as an email.
    """

    logger.info(f"{ont_user_run.run_name}: Touching StatusDB...")
    ont_user_run.touch_db_entry()

    if not ont_user_run.is_synced():
        logger.info(f"{ont_user_run.run_name}: Run is not fully synced, skipping.")
    else:

        if ont_user_run.is_transferred():
            logger.warning(
                f"{ont_user_run.run_name}: Run is already logged as transferred, skipping."
            )
        else:
            logger.info(f"{ont_user_run.run_name}: Processing transfer...")

            # Assert all files are in place
            logger.info(f"{ont_user_run.run_name}: Asserting run contents...")
            ont_user_run.assert_contents()

            # Update StatusDB
            logger.info(f"{ont_user_run.run_name}: Updating StatusDB...")
            ont_user_run.update_db_entry()

            # Copy HTML report
            logger.info(f"{ont_user_run.run_name}: Putting HTML report on GenStat...")
            ont_user_run.copy_html_report()

            # Copy metadata
            logger.info(f"{ont_user_run.run_name}: Copying metadata...")
            ont_user_run.copy_metadata()

            # Transfer run
            logger.info(f"{ont_user_run.run_name}: Transferring to cluster...")
            ont_user_run.transfer_run()

            # Update transfer log
            logger.info(f"{ont_user_run.run_name}: Updating transfer log...")
            ont_user_run.update_transfer_log()

            # Archive run
            logger.info(f"{ont_user_run.run_name}: Archiving run...")
            ont_user_run.archive_run()


def process_qc_run(ont_qc_run: ONT_qc_run):
    """This control function orchestrates the sequential execution of the ONT_qc_run class methods.

    For a single ONT QC run...

        - Ensure there is a database entry corresponding to an ongoing run

        If not fully synced:
            - Skip

        If fully synced:
            - Ensure all necessary files to proceed with processing are present
            - Update the StatusDB entry
            - Copy HTML report to GenStat

            If Anglerfish has not been run:
                If Anglerfish is ongoing:
                    - Skip

                If Anglerfish is not ongoing:
                    - Run Anglerfish

            If Anglerfish has been run:
                If Anglerfish failed:
                    - Throw error

                If Anglerfish finished successfully:
                    - Copy metadata
                    - Transfer run to cluster
                    - Update transfer log
                    - Archive run

    Any errors raised here-in should be sent with traceback as an email.
    """

    logger.info(f"{ont_qc_run.run_name}: Touching StatusDB...")
    ont_qc_run.touch_db_entry()

    # Is the run fully synced?
    if not ont_qc_run.is_synced():
        logger.info(f"{ont_qc_run.run_name}: Run is not fully synced, skipping.")
    else:

        # Assert all files are in place
        logger.info(f"{ont_qc_run.run_name}: Asserting run contents...")
        ont_qc_run.assert_contents()

        # Update StatusDB
        logger.info(f"{ont_qc_run.run_name}: Updating StatusDB...")
        ont_qc_run.update_db_entry()

        # Copy HTML report
        logger.info(f"{ont_qc_run.run_name}: Putting HTML report on GenStat...")
        ont_qc_run.copy_html_report()

        # Anglerfish
        logger.info(
            f"{ont_qc_run.run_name}: Checking whether Anglerfish has been run..."
        )

        anglerfish_exit_code = ont_qc_run.get_anglerfish_exit_code()

        if isinstance(anglerfish_exit_code, int) and anglerfish_exit_code > 0:
            logger.warning(
                f"{ont_qc_run.run_name}: Anglerfish has failed, throwing error."
            )
            raise AssertionError(f"{ont_qc_run.run_name}: Anglerfish failed.")

        elif anglerfish_exit_code is None:
            logger.info(
                f"{ont_qc_run.run_name}: Anglerfish has not been run, continuing."
            )

            logger.info(
                f"{ont_qc_run.run_name}: Checking whether Anglerfish is ongoing..."
            )

            anglerfish_pid = ont_qc_run.get_anglerfish_pid()
            if anglerfish_pid:
                logger.info(
                    f"{ont_qc_run.run_name}: Anglerfish is ongoing with process ID {anglerfish_pid}, skipping."
                )
            else:
                logger.info(
                    f"{ont_qc_run.run_name}: Anglerfish is not ongoing, continuing."
                )

                logger.info(
                    f"{ont_qc_run.run_name}: Fetching Anglerfish samplesheet..."
                )
                if not ont_qc_run.fetch_anglerfish_samplesheet():
                    logger.info(
                        f"{ont_qc_run.run_name}: Could not find Anglerfish sample sheet, skipping."
                    )
                else:
                    logger.info(f"{ont_qc_run.run_name}: Starting Anglerfish...")
                    ont_qc_run.run_anglerfish()

        elif isinstance(anglerfish_exit_code, int) and anglerfish_exit_code == 0:
            logger.info(
                f"{ont_qc_run.run_name}: Anglerfish has finished successfully, continuing."
            )

            if ont_qc_run.is_transferred():
                logger.warning(
                    f"{ont_qc_run.run_name}: Run is already logged as transferred, skipping."
                )

            else:
                logger.info(f"{ont_qc_run.run_name}: Processing transfer...")

                # Copy metadata
                logger.info(f"{ont_qc_run.run_name}: Copying metadata...")
                ont_qc_run.copy_metadata()

                # Transfer run
                logger.info(f"{ont_qc_run.run_name}: Transferring to cluster...")
                ont_qc_run.transfer_run()

                # Update transfer log
                logger.info(f"{ont_qc_run.run_name}: Updating transfer log...")
                ont_qc_run.update_transfer_log()

                # Archive run
                logger.info(f"{ont_qc_run.run_name}: Archiving run...")
                ont_qc_run.archive_run()


def process_run(run_abspath: str):
    """This gate function instantiates an appropriate subclass
    for the ONT run and processes it"""

    ont_run = ONT_run(run_abspath)

    if ont_run.run_type == "user_run":
        process_user_run(ONT_user_run(run_abspath))
    elif ont_run.run_type == "qc_run":
        process_qc_run(ONT_qc_run(run_abspath))
    else:
        raise AssertionError("Run type invalid.")


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

            data_dirs = CONFIG["nanopore_analysis"]["run_types"][run_type]["data_dirs"]
            ignore_dirs = CONFIG["nanopore_analysis"]["run_types"][run_type][
                "ignore_dirs"
            ]

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
    ont_run.update_db_entry(force_update=True)
