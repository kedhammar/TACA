"""Analysis methods for sequencing runs produced by Element instruments."""

import glob
import logging
import os

from taca.element.Element_Runs import Aviti_Run
from taca.utils.config import CONFIG

logger = logging.getLogger(__name__)


def _upload_to_statusdb(run):
    """Triggers the upload to statusdb.

    :param Run run: the object run
    """
    pass


def run_preprocessing(given_run):
    """Run demultiplexing in all data directories.

    :param str given_run: Process a particular run instead of looking for runs
    """

    def _process(run):
        """Process a run/flowcell and transfer to analysis server.

        :param taca.element.Run run: Run to be processed and transferred
        """
        # TODO: Fetch statusdb document for run
        # TODO: Get previous status of run from statusdb document
        sequencing_done = run.check_sequencing_status()
        demultiplexing_status = run.get_demultiplexing_status()
        if not sequencing_done:
            # TODO: compare previous status with current status and update statusdb document if different
            return
        elif sequencing_done and demultiplexing_status == "not started":
            if not run.manifest_exists():  # Assumes that we use the same manifest as for sequencing. TODO: demux settings need to be added to the original manifest by lims
                # TODO: email operator that manifest is missing
                return
            # Start demux
            run.start_demux()
            # TODO: compare previous status with current status and update statusdb document if different
            return
        elif sequencing_done and demultiplexing_status == "ongoing":
            # TODO: compare previous status with current status and update statusdb document if different
            return
        elif sequencing_done and demultiplexing_status == "finished":
            # Sync metadata to ngi-data-ns
            # check if run is transferred or transfer is ongoing
            # if run has not been transferred and transfer is not ongoing
            # make a hidden file to indicate that transfer has started
            # compare previous status with current status and update statusdb document if different
            # Also update statusdb with a timestamp of when the transfer started
            # transfer run to miarka
            # remove hidden file if transfer was successful
            # Update transfer log
            # update statusdb document
            # archive run to nosync
            # update statusdb document
            # elif run is being transferred (hidden file exists)
            # compare previous status with current status and update statusdb document if different
            # return
            # elif run is already transferred (in transfer log)
            # compare previous status with current status and update statusdb document if different
            # warn that transferred run has not been archived
            pass

    if given_run:
        run = Aviti_Run(
            run
        )  # TODO: Needs to change if more Element machines are aquired in the future
        _process(runObj)
    else:
        data_dirs = CONFIG.get("element_analysis").get(
            "data_dirs"
        )  # TODO: add to config
        for data_dir in data_dirs:
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
                    pass
