"""
Nanopore analysis methods for TACA
"""
import os
import logging
import glob
import csv
import subprocess
import shutil
import smtplib

from datetime import datetime
from taca.utils.config import CONFIG
from taca.utils.transfer import RsyncAgent
from taca.utils.misc import send_mail

logger = logging.getLogger(__name__)

def find_runs_to_process():
    data_dirs = CONFIG.get('nanopore_analysis').get('data_dirs')
    found_run_dirs = []
    for data_dir in data_dirs:
        try:
            found = [os.path.join(data_dir, run_dir) for run_dir in os.listdir(data_dir)
                     if os.path.isdir(os.path.join(data_dir, run_dir))
                     and run_dir != 'nosync']
        except OSError:
            logger.warn("There was an issue locating the following directory: " + data_dir +
                        ". Please check that it exists and try again.")
        if found:
            for found_dir in found:
                found_run_dirs.append(found_dir)
        else:
            logger.warn("Could not find any run directories in " + data_dir + ". Skipping.")
    return found_run_dirs

def process_run(run_dir):
    logger.info("Processing run: " + run_dir)
    summary_file = os.path.join(run_dir, "final_summary.txt")
    demux_dir = os.path.join(run_dir, "nanoseq_output")  # TODO: Check actual name of output dir
    sample_sheet = os.path.join(run_dir, "sample_sheet.csv")  # TODO: Check actual name for sample sheet
    analysis_exit_status_file = os.path.join(run_dir, ".exitcode_for_nanoseq")
    email_recipients = CONFIG.get('mail').get('recipients')
    if os.path.isfile(summary_file) and not os.path.isdir(demux_dir):
        logger.info("Sequencing done for run " + run_dir + ". Attempting to start analysis.")
        if os.path.isfile(sample_sheet):
            start_analysis_pipeline(run_dir)
#        else if data available in lims:
#            make sample_sheet
#            start_analysis_pipeline(run_dir)
        else:
            logger.warn("Samplesheet not found for run " + run_dir + ". Operator notified. Skipping.")
            email_subject = ("Samplesheet missing for run {}".format(os.path.basename(run_dir)))
            email_message = """The samplesheet for run {run} is missing and the
            information can't be found in LIMS. Please add the samplesheet to {run}.""".format(run=run_dir)
            send_mail(email_subject, email_message, email_recipients)
    elif os.path.isdir(demux_dir) and not os.path.isfile(analysis_exit_status_file):
        logger.info("Analysis has started for run " + run_dir +" but is not yet done. Skipping.")
    elif os.path.isdir(demux_dir) and os.path.isfile(analysis_exit_status_file):
        analysis_successful = check_exit_status(analysis_exit_status_file)
        if analysis_successful:
            run_id = os.path.basename(run_dir)
            transfer_file_name = CONFIG.get('nanopore_analysis').get('transfer').get('transfer_file')
            transfer_log = os.path.join(os.path.dirname(run_dir), transfer_file_name)
            if is_not_transferred(run_id, transfer_log):
                transfer_run(run_dir)
                update_transfer_log(run_id, transfer_log)
                logger.info("Run "+ run_dir + " has been synced to the analysis cluster.")
                archive_run(run_dir)
                logger.info("Run " + run_dir + " is finished and has been archived.")
                email_subject = ("Run successfully processed: {}".format(os.path.basename(run_dir)))
                email_message = """Run {} has been analysed, transferred and archived
                successfully.""".format(run_dir)
                send_mail(email_subject, email_message, email_recipients)
        else:
            logger.warn("Analysis pipeline exited with a non-zero exit status for run " + run_dir + ". Notifying operator.")
            email_subject = ("Analysis failed for run {}".format(os.path.basename(run_dir)))
            email_message = """The analysis failed for run {run}.
            Please review the logfiles in {run}.""".format(run=run_dir)
            send_mail(email_subject, email_message, email_recipients)
    else:
        logger.info("Run " + run_dir + " not finished yet. Skipping.")
    return

def start_analysis_pipeline(run_dir):
    # start analysis detatched
    analysis_command = "nextflow run nf-core/nanoseq -r dev --help ; echo $? > .exitcode_for_taca.txt"  # TODO: Change to actual command
    try:
        p_handle = subprocess.Popen(analysis_command, stdout=subprocess.PIPE, shell=True, cwd=run_dir)
        logger.info("Started analysis for run " + run_dir)
    except subprocess.CalledProcessError:
        logger.warn("An error occurred while starting the analysis for run " + run_dir + ". Please check the logfile for info.")
    return

def check_exit_status(status_file):
    # Read pipeline exit status file and return True if 0, False if anything else
    with open(status_file, 'r') as f:
        exit_status = f.readline().strip()
    if exit_status == '0':
        return True
    else:
        return False

def is_not_transferred(run_id, transfer_log):
    # Return True if run id not in transfer.tsv, else False
        with open(transfer_log, 'r') as f:
            if run_id not in f.read():
                return True
            else:
                return False

def transfer_run(run_dir):
    #rsync dir to irma
    logger.info("Transferring run " + run_dir + " to analysis cluster")
    destination = CONFIG.get("nanopore_analysis").get("transfer").get("destination")
    rsync_opts = {"--no-o" : None, "--no-g" : None, "--chmod" : "g+rw", "-r" : None}
    connection_details = CONFIG.get("nanopore_analysis").get("transfer").get("analysis_server")
    transfer_object = RsyncAgent(run_dir, dest_path=destination, remote_host=connection_details["host"], remote_user=connection_details["user"], validate=False, opts=rsync_opts)
    try:
        transfer_object.transfer()
    except RsyncError:
        logger.warn("An error occurred while transferring " + run_dir + " to the ananlysis server. Please check the logfiles")
    return

def update_transfer_log(run_id, transfer_log):
    try:
        with open(transfer_log, 'a') as f:
            tsv_writer = csv.writer(f, delimiter='\t')
            tsv_writer.writerow([run_id, str(datetime.now())])
    except IOError:
        logger.warn("Could not update the transfer logfile for run " + run_id + ". Please make sure " + transfer_log + " gets updated.")
    return

def archive_run(run_dir):
    # mv dir to nosync
    logger.info("Archiving run " + run_dir)
    archive_dir = CONFIG.get("nanopore_analysis").get("finished_dir")
    archive_path = os.path.join(os.path.dirname(run_dir), archive_dir, os.path.basename(run_dir))
    try:
        shutil.move(run_dir, archive_path)
        logger.info("Successfully archived " + run_dir)
    except shutil.Error:
        logger.warn("An error occurred when archiving " + run_dir + ". Please check the logfile for more info.")
    return

def run_preprocessing(run):
    if run:
        process_run(run)
    else:
        runs_to_process = find_runs_to_process()
        for run_dir in runs_to_process:
            process_run(run_dir)
