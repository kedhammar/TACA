"""Analysis methods for sequencing runs produced by Element instruments."""

import glob
import logging
import os

from taca.element.Aviti_Runs import Aviti_Run
from taca.utils.config import CONFIG

logger = logging.getLogger(__name__)


def run_preprocessing(given_run):
    """Run demultiplexing in all data directories.

    :param str given_run: Process a particular run instead of looking for runs
    """

    def _process(run):
        """Process a run/flowcell and transfer to analysis server.

        :param taca.element.Run run: Run to be processed and transferred
        """
        try:
            run.parse_run_parameters()
        except FileNotFoundError:
            logger.warn(
                f"Cannot reliably set NGI_run_id for {run} due to missing RunParameters.json. Aborting run processing"
            )
            raise

        sequencing_done = run.check_sequencing_status()
        demultiplexing_status = run.get_demultiplexing_status()
        if not sequencing_done:  # Sequencing ongoing
            run.status = "sequencing"
            if run.status_changed:
                run.update_statusdb()
        elif (
            sequencing_done and demultiplexing_status == "not started"
        ):  # Sequencing done. Start demux
            if (
                not run.manifest_exists()
            ):  # TODO: this should check for the zip file in lims output location
                logger.warn(
                    f"Run manifest is missing for {run}, demultiplexing aborted"
                )
                # TODO: email operator warning
                return
            elif run.manifest_exists():
                os.mkdir(run.demux_dir)
                run.copy_manifests()
                run_manifests = glob.glob(
                    os.path.join(
                        run.run_dir, "RunManifest_*.csv"
                    )  # TODO: is this filename right?
                )
                sub_demux_count = 0
                for run_manifest in run_manifests.sort():
                    if len(run_manifests) == 1:
                        demux_dir = run.demux_dir
                    elif len(run_manifests) > 1:
                        demux_dir = f"Demultiplexing_{sub_demux_count}"
                    os.mkdir(demux_dir)
                    run.start_demux(run_manifest, demux_dir)
                    sub_demux_count += 1
                run.status = "demultiplexing"
                if run.status_changed:
                    run.update_statusdb()
        elif sequencing_done and demultiplexing_status == "ongoing":
            run.status = "demultiplexing"
            if run.status_changed:
                run.update_statusdb()
            return
        elif sequencing_done and demultiplexing_status == "finished":
            transfer_status = run.get_transfer_status()
            if transfer_status == "not started":
                run.aggregate_demux_results() # TODO: if multiple demux dirs, aggregate the results into Demultiplexing?
                run.sync_metadata()
                run.make_transfer_indicator()
                run.status = "transferring"
                if run.status_changed:
                    run.update_statusdb()
                    # TODO: Also update statusdb with a timestamp of when the transfer started
                run.transfer()  # I think this should be a detached command as well
            elif transfer_status == "ongoing":
                run.status = "transferring"
                if run.status_changed:
                    run.update_statusdb()
                logger.info(f"{run} is being transferred. Skipping.")
                return
            elif transfer_status == "rsync done":
                if run.rsync_success():
                    run.remove_transfer_indicator()
                    run.update_transfer_log()
                    run.status = "transferred"
                    if run.status_changed:
                        run.update_statusdb()
                    run.archive()
                    run.status = "archived"
                    if run.status_changed:
                        run.update_statusdb()
                else:
                    run.status = "transfer failed"
                    logger.warning(f"An issue occurred while transfering {run} to the analysis cluster." )
                    # TODO: email warning to operator
            elif transfer_status == "unknown":
                logger.warning(
                    f"The run {run} has already been transferred but has not been archived. Please investigate"
                )
                # TODO: email operator warning
                return
            else:
                logger.warning(f"Unknown transfer status of run {run}. Please investigate")

    if given_run:
        run = Aviti_Run(given_run)
        # TODO: Needs to change if more types of Element machines are aquired in the future

        _process(run)
    else:
        data_dirs = CONFIG.get("element_analysis").get(
            "data_dirs"
        )  # TODO: add to config
        for data_dir in data_dirs:  # TODO: make sure to look in both side A and B
            # Run folder looks like DATE_*_*_*, the last section is the FC name.
            runs = glob.glob(
                os.path.join(data_dir, "[1-9]*_*_*_*")
            )  # TODO: adapt to aviti format
            for run in runs:
                runObj = Aviti_Run(run)
                try:
                    _process(runObj)
                except:  # TODO: chatch error message and print it
                    # This function might throw and exception,
                    # it is better to continue processing other runs
                    logger.warning(f"There was an error processing the run {run}")
                    # TODO: Think about how to avoid silent errors (email?)
                    pass
