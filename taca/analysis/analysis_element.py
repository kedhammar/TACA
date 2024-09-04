"""Analysis methods for sequencing runs produced by Element instruments."""

import glob
import logging
import os

from taca.element.Element_Runs import Aviti_Run
from taca.utils.config import CONFIG
from taca.utils import statusdb


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
        #TODO: Fetch statusdb document for run
        #TODO: Get previous status of run from statusdb document
        sequencing_done = run.check_sequencing_status()
        demultiplexing_status = run.get_demultiplexing_status()
        if not sequencing_done:
            #TODO: compare previous status with current status and update statusdb document if different
            return
        elif sequencing_done and demultiplexing_status == "not started":
            if not run.manifest_exists():
                #TODO: email operator warning
                return
            elif run.manifest_exists():
                # Get sample info from manifest
                sample_info = run.get_sample_info()
                sample_types = run.get_sample_types(sample_info)
                if len(sample_types) == 1:
                    run.start_demux()
                elif len(sample_types) > 1:
                    for sample_type in sample_types:
                        run.make_manifest(sample_info, sample_type)
                        run.start_demux()
                else:
                    #TODO: warn that no samples were found in the run manifest
                    return
                #TODO: compare previous status with current status and update statusdb document if different
        elif sequencing_done and demultiplexing_status == "ongoing":
            #TODO: compare previous status with current status and update statusdb document if different
            return
        elif sequencing_done and demultiplexing_status == "finished":
            transfer_file = CONFIG.get('Element').get('Aviti').get('transfer_log')
            if not run.is_transferred(transfer_file) and not run.transfer_ongoing():
                run.sync_metadata()
                run.make_transfer_indicator()
                #TODO: compare previous status with current status and update statusdb document if different
                    # Also update statusdb with a timestamp of when the transfer started
                run.transfer()
                run.remove_transfer_indicator()
                run.update_transfer_log(transfer_file)
                #TODO: update statusdb document
                run.archive()
            elif not run.is_transferred(transfer_file) and run.transfer_ongoing():
                #TODO: compare previous status with current status and update statusdb document if different
                logger.info("Run is being transferred. Skipping.")
                return
            elif run.is_transferred(transfer_file):
                #TODO: compare previous status with current status and update statusdb document if different
                # warn that transferred run has not been archived
                logger.warn("The run has already been transferred but has not been archived. Please investigate")
                return
            else:
                logger.warn("Unknown transfer status. Please investigate")
                
            

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
