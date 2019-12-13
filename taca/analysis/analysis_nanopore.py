"""
Nanopore analysis methods for TACA
"""
import os
import logging
import glob
from taca.utils.config import CONFIG

logger = logging.getLogger(__name__)

def find_runs_to_process():
    data_dirs = CONFIG.get('nanopore_analysis').get('data_dirs')
    found_run_dirs = []
    for data_dir in data_dirs:
        found = [os.path.join(data_dir, run_dir) for run_dir in os.listdir(data_dir)
                              if os.path.isdir(os.path.join(data_dir, run_dir))]
        for found_dir in found:
            found_run_dirs.append(found_dir)  # TODO: See if there is a neater way to do this
    return found_run_dirs

def process_run(run_dir):
    logger.info("Processing run: " + run_dir)
    run_id = os.path.basename(run_dir)
    summary_file = os.path.join(run_dir, "final_summary.txt")
    demux_dir = os.path.join(run_dir, "nanoseq_output")  # TODO: Check actual name of output dir
    sample_sheet = os.path.join(run_dir, "sample_sheet.csv")  # TODO: Check actual name for sample sheet
    analysis_exit_status_file = os.path.join(run_dir, ".exitcode_for_nanoseq")
    if os.path.isfile(summary_file) and not os.path.isdir(demux_dir):
        logger.info("Sequencing done for run " + run_dir + ". Attempting to start analysis.")
        if os.path.isfile(sample_sheet):
            start_analysis_pipeline(run_dir)
#        else if data available in lims:
#            make sample_sheet
#            start_analysis_pipeline(run_dir)
        else:  # Samplesheet missing and info not found in lims. Notify user and skip.
            logger.warn("Samplesheet not found for run " + run_dir + ". Operator notified. Skipping.")
    elif os.path.isdir(demux_dir) and not os.path.isfile(analysis_exit_status_file):
        logger.info("Analysis has started for run " + run_dir +" but is not yet done. Skipping.")
    elif os.path.isdir(demux_dir) and os.path.isfile(analysis_exit_status_file):
        analysis_successful = check_exit_status(analysis_exit_status_file)
        if analysis_successful:
            if is_not_transferred(run_id):
                transfer_run(run_dir)
            else:
                archive_run(run_dir)
                logger.info("Run " + run_dir + " is finished and has been archived.")
                # TODO: email operator
        else:
            logger.warn("Analysis pipeline exited with a non-zero exit status for run " + run_dir + ". Notifying operator.")
            # TODO: Email operator
    else:
        logger.info("Run " + run_dir + " not finished yet. Skipping.")
    return

def start_analysis_pipeline(run_dir):
    # start analysis detatched
    logger.info("Started analysis for run " + run_dir)
    return

def check_exit_status(status_file):
    #read file and return True if 0, False if anything else
    return True

def is_not_transferred(run_id):
    #return True if id not in transfer.tsv, else False
    return True

def transfer_run(run_dir):
    #rsync dir to irma
    logger.info("Transferring run " + run_dir + " to analysis cluster")
    return

def archive_run(run_dir):
    # mv dir to nosync
    logger.info("Archiving run " + run_dir)
    return

def run_preprocessing(run):
    if run:
        process_run(run)
    else:
        runs_to_process = find_runs_to_process()
        for run_dir in runs_to_process:
            process_run(run_dir)
