"""Analysis methods for sequencing runs produced by Element instruments."""

import glob
import logging
import os

from taca.element.Element_Runs import Aviti_Run
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
        # Check if sequencing is finished. (is the final file there and was it completed OK)
        # if sequencing is not done
            # Update statusdb?
            # return
        # else If sequencing finished and demux not started
            # Update statusdb
            # Get/generate sample sheet
            # Start demux
        # else if sequencing finished and demux ongoing
            # do nothing
        # Else if sequencing started and demux finished
            # check if run is transferred or transfer is ongoing
            # if run has not been transferred and transfer is not ongoing
                # make a hidden file to indicate that transfer has started
                # transfer run to miarka 
                # remove hidden file if transfer was successful
                # Update transfer log
                # archive run to nosync
            # elif run is being transferred (hidden file exists)
                # return
            # elif run is already transferred (in transfer log)
                # warn that transferred run has not been archived

        

    if given_run:
        run = Aviti_Run(run) #TODO: Needs to change if more Element machines are aquired in the future
        _process(runObj)
    else:
        data_dirs = CONFIG.get("element_analysis").get("data_dirs") #TODO: add to config
        for data_dir in data_dirs:
            # Run folder looks like DATE_*_*_*, the last section is the FC name.
            runs = glob.glob(os.path.join(data_dir, "[1-9]*_*_*_*")) #TODO: adapt to aviti format
            for run in runs:
                runObj = Aviti_Run(run)
                try:
                    _process(runObj)
                except: #TODO: chatch error message and print it
                    # This function might throw and exception,
                    # it is better to continue processing other runs
                    logger.warning(f"There was an error processing the run {run}")
                    pass
