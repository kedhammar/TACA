"""Nanopore analysis methods for TACA."""
import os
import logging

from dateutil.parser import parse
from taca.utils.config import CONFIG
from taca.utils.misc import send_mail
from taca.nanopore.ont_transfer_class import PromethionTransfer, MinionTransfer

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


def ont_updatedb_from_cli(run_abspath):

    if is_date(os.path.basename(run_abspath).split("_")[0]):
        if "minion" in run_abspath:
            ont_run = MinionTransfer(os.path.abspath(run_abspath))
        elif "promethion" in run_abspath:
            ont_run = PromethionTransfer(os.path.abspath(run_abspath))
        ont_run.update_db(force_update=True)
    else:
        logger.warning(
            "The specified path is not a flow cell. Please "
            "provide the full path to the flow cell you wish to process."
        )
