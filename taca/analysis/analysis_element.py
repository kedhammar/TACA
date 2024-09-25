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
            logger.warning(
                f"Cannot reliably set NGI_run_id for {run} due to missing RunParameters.json. Aborting run processing"
            )
            raise

        #### Sequencing status ####
        sequencing_done = run.check_sequencing_status()
        if not sequencing_done:  # Sequencing ongoing
            run.status = "sequencing"
            if run.status_changed:
                run.update_statusdb()
            return

        #### Demultiplexing status ####
        demultiplexing_status = run.get_demultiplexing_status()
        if demultiplexing_status == "not started":
            # Sequencing done. Start demux
            if run.manifest_exists():
                os.mkdir(run.demux_dir)
                run.copy_manifests()
                run_manifests = glob.glob(
                    os.path.join(
                        run.run_dir, "RunManifest_*.csv"
                    )  # TODO: is this filename right?
                )
                sub_demux_count = 0
                for run_manifest in run_manifests.sort():
                    demux_dir = f"Demultiplexing_{sub_demux_count}"
                    os.mkdir(demux_dir)
                    run.start_demux(run_manifest, demux_dir)
                    sub_demux_count += 1
                run.status = "demultiplexing"
                if run.status_changed:
                    run.update_statusdb()
                return
            else:
                logger.warning(
                    f"Run manifest is missing for {run}, demultiplexing aborted"
                )
                # TODO: email operator warning
                return
        elif demultiplexing_status == "ongoing":
            run.status = "demultiplexing"
            if run.status_changed:
                run.update_statusdb()
            return
          
        elif demultiplexing_status != "finished":
            logger.warning(
                f"Unknown demultiplexing status {demultiplexing_status} of run {run}. Please investigate"
            )
            return

        #### Transfer status ####
        transfer_status = run.get_transfer_status()
        if transfer_status == "not started":
            demux_results_dirs = glob.glob(
                os.path.join(run.run_dir, "Delmultiplexing_*")
            )
            run.aggregate_demux_results(demux_results_dirs)
            run.upload_demux_results_to_statusdb()
            run.sync_metadata()
            run.make_transfer_indicator()
            run.status = "transferring"
            if run.status_changed:
                run.update_statusdb()
                # TODO: Also update statusdb with a timestamp of when the transfer started
            run.transfer()
            return
        elif transfer_status == "ongoing":
            run.status = "transferring"
            if run.status_changed:
                run.update_statusdb()
            logger.info(f"{run} is being transferred. Skipping.") # TODO: fix formatting, currently prints "ElementRun(20240910_AV242106_B2403418431) is being transferred"
            return
        elif transfer_status == "rsync done":
            if run.rsync_successful():
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
                logger.warning(
                    f"An issue occurred while transfering {run} to the analysis cluster."
                )
                # TODO: email warning to operator
            return
        elif transfer_status == "unknown":
            logger.warning(
                f"The run {run} has already been transferred but has not been archived. Please investigate"
            )
            # TODO: email operator warning
            return
        else:
            # TODO Merge with the one above?
            logger.warning(
                f"Unknown transfer status {transfer_status} of run {run}. Please investigate"
            )
            return

    if given_run:
        run = Aviti_Run(given_run, CONFIG)
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
                runObj = Aviti_Run(run, CONFIG)
                try:
                    _process(runObj)
                except:  # TODO: chatch error message and print it
                    # This function might throw and exception,
                    # it is better to continue processing other runs
                    logger.warning(f"There was an error processing the run {run}")
                    # TODO: Think about how to avoid silent errors (email?)
                    pass
