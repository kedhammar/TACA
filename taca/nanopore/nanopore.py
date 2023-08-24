import os
import logging
import csv
import shutil
import glob
import pathlib
import re
import json
import pandas as pd

from taca.utils.statusdb import NanoporeRunsConnection
from datetime import datetime
from taca.utils.config import CONFIG
from taca.utils.transfer import RsyncAgent, RsyncError

logger = logging.getLogger(__name__)

class Nanopore(object):
    """General Nanopore run"""
    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.run_id = os.path.basename(run_dir)
        self.summary_file = glob.glob(run_dir + '/final_summary*.txt')
        self.top_dir = str(pathlib.Path(self.run_dir).parent.parent)
        self.experiment_id = os.path.basename(self.top_dir)

    def is_not_transferred(self):
        """Return True if run id not in transfer.tsv, else False."""
        with open(self.transfer_log, 'r') as f:
            return self.run_id not in f.read()

    def transfer_run(self):
        """rsync dir to destination specified in config file."""
        destination = self.transfer_details.get('destination')
        rsync_opts = self.transfer_details.get('rsync_options')
        for k, v in rsync_opts.items():
            if v == 'None':
                rsync_opts[k] = None
        connection_details = self.transfer_details.get('analysis_server', None)
        logger.info('Transferring run {} to {}'.format(self.run_id, 
                                                       connection_details['host'] if connection_details 
                                                       else destination))
        if connection_details:
            transfer_object = RsyncAgent(self.run_dir,
                                        dest_path=destination,
                                        remote_host=connection_details.get('host'),
                                        remote_user=connection_details.get('user'),
                                        validate=False,
                                        opts=rsync_opts)
        else:
            transfer_object = RsyncAgent(self.run_dir,
                                        dest_path=destination,
                                        validate=False,
                                        opts=rsync_opts)
        try:
            transfer_object.transfer()
        except RsyncError:
            logger.warn('An error occurred while transferring {} to the '
                        'analysis server. Please check the logfiles'.format(self.run_dir))
            return False
        return True

    def transfer_metadata(self):
        """Copies run dir (excluding seq data) to metadata dir"""

        exclude_patterns = [
            # Main seq dirs
            "**/bam*/***",
            "**/fast5*/***",
            "**/fastq*/***",
            # Any files found elsewhere
            "*.bam*",
            "*.fast5*",
            "*.fastq*",
        ]

        exclude_patterns_quoted = ["'" + pattern + "'" for pattern in exclude_patterns]

        src = self.run_dir
        dst = os.path.join(self.metadata_dir)

        os.system(
            f"rsync -rv --exclude={{{','.join(exclude_patterns_quoted)}}} {src} {dst}"
        )

    def transfer_html_report(self):

        glob_html = glob.glob(self.run_dir + "/report*.html")

        if len(glob_html) == 0:
            error_message = f"Run {self.run_id} is marked as finished, but missing .html report file."
            logger.error(error_message)
            raise AssertionError(error_message)
        elif len(glob_html) > 1:
            error_message = f"Run {self.run_id} is marked as finished, but contains conflicting .html report files."
            logger.error(error_message)
            raise AssertionError(error_message)

        logger.info(f"Transferring the run report to ngi-internal.")

        # Transfer the MinKNOW .html report file to ngi-internal, renaming it to the full run ID. Requires password-free SSH access.
        report_dest_path = os.path.join(
            CONFIG["nanopore_analysis"]["ont_transfer"]["minknow_reports_dir"],
            f"report_{self.run_id}.html",
        )
        transfer_object = RsyncAgent(
            glob_html[0],
            dest_path=report_dest_path,
            validate=False,
        )
        try:
            transfer_object.transfer()
            logger.info(
                f"Successfully transferred the MinKNOW report of run {self.run_id}"
            )
        except RsyncError:
            msg = f"An error occurred while attempting to transfer the report {glob_html[0]} to {report_dest_path}"
            logger.error(msg)
            raise RsyncError(msg)

    def update_db(self, force_update=False):
        """Check run vs statusdb. Create or update run entry."""

        logger.info("Updating database with run {}".format(self.run_id))

        run_pattern = re.compile(
            "^(\d{8})_(\d{4})_([0-9a-zA-Z]+)_([0-9a-zA-Z]+)_([0-9a-zA-Z]+)$"
        )

        sesh = NanoporeRunsConnection(CONFIG["statusdb"], dbname="nanopore_runs")

        if re.match(run_pattern, self.run_id):
            logger.info(f"Run {self.run_id} looks like a run directory, continuing.")
        else:
            error_message = f"Run {self.run_id} does not match the regex of a run directory (yyyymmdd_hhmm_pos|device_fcID_hash)."
            logger.error(error_message)
            raise AssertionError(error_message)

        # If no run document exists in the database, ceate an ongoing run document
        if not sesh.check_run_exists(self):
            logger.info(
                f"Run {self.run_id} does not exist in the database, creating entry for ongoing run."
            )

            run_path_file = os.path.join(self.run_dir, "run_path.txt")
            assert os.path.isfile(run_path_file), f"Couldn't find {run_path_file}"

            pore_count_history_file = os.path.join(
                self.run_dir, "pore_count_history.csv"
            )
            assert os.path.isfile(
                pore_count_history_file
            ), f"Couldn't find {pore_count_history_file}"

            sesh.create_ongoing_run(self, run_path_file, pore_count_history_file)
            logger.info(f"Successfully created db entry for ongoing run {self.run_id}.")

        # If the run document is marked as "ongoing" or database is being manually updated
        if sesh.check_run_status(self) == "ongoing" or force_update == True:
            logger.info(
                f"Run {self.run_id} exists in the database with run status: {sesh.check_run_status(self)}"
            )

            # If the run is finished
            if len(self.summary_file) != 0:

                logger.info(
                    f"Run {self.run_id} has finished sequencing, updating the db entry."
                )

                # Instantiate json (dict) to update the db with
                db_update = {}

                # Parse report_*.json
                try:
                    self.parse_minknow_json(db_update)
                except BaseException as e:
                    logger.error(f"Failed parse_minknow_json() for run {self.run_id}")
                    raise e

                # Parse pore_activity_*.csv
                try:
                    self.parse_pore_activity(db_update)
                except BaseException as e:
                    logger.error(f"Failed parse_pore_activity() for run {self.run_id}")
                    raise e

                # Update the DB entry
                try:
                    sesh.finish_ongoing_run(self, db_update)
                    logger.info(
                        f"Successfully updated the db entry of run {self.run_id}"
                    )
                except BaseException as e:
                    logger.error(f"Failed finish_ongoing_run() for run {self.run_id}")
                    raise e

                # Transfer the MinKNOW run report
                try:
                    self.transfer_html_report()
                except BaseException as e:
                    logger.error(f"Failed transfer_html_report() for run {self.run_id}")
                    raise e

            else:
                logger.info(
                    f"Run {self.run_id} has not finished sequencing, do nothing."
                )

        # if the run document is marked as "finished"
        elif sesh.check_run_status(self) == "finished":
            logger.info(
                f"Run {self.run_id} exists in the database as an finished run, do nothing."
            )

    def parse_pore_activity(self, db_update):

        logger.info(f"Parsing pore activity of run {self.run_id}")

        pore_activity = {}

        # Find the pore activity file
        pore_activity_filename = "pore_activity_*.csv"

        glob_pore_activity = glob.glob(self.run_dir + f"/{pore_activity_filename}")

        if len(glob_pore_activity) == 0:
            error_message = f"Run {self.run_id} is marked as finished, but missing {pore_activity_filename}"
            logger.error(error_message)
            raise AssertionError(error_message)
        elif len(glob_pore_activity) > 1:
            error_message = f"Run {self.run_id} is marked as finished, but contains conflicting {pore_activity_filename} files."
            logger.error(error_message)
            raise AssertionError(error_message)

        # Use pandas to pivot the data into a more manipulable dataframe
        df = pd.read_csv(glob_pore_activity[0])
        df.sort_values(by="Experiment Time (minutes)", inplace=True)
        df = df.pivot_table(
            "State Time (samples)", "Experiment Time (minutes)", "Channel State"
        )

        # Use pore counts to calculate new metrics
        df["all"] = df.sum(axis=1)
        df["healthy"] = df.strand + df.adapter + df.pore
        df["productive"] = df.strand + df.adapter

        df["health"] = df["healthy"] / df["all"]
        df["efficacy"] = df["productive"] / df["healthy"]

        # Look at peaks within 1st hour of the run and define some metrics
        df_h1 = df[0:60]
        pore_activity["peak_pore_health_pc"] = round(
            100 * float(df_h1.loc[df_h1.health == df_h1.health.max(), "health"]), 2
        )
        pore_activity["peak_pore_efficacy_pc"] = round(
            100 * float(df_h1.loc[df_h1.efficacy == df_h1.efficacy.max(), "efficacy"]),
            2,
        )

        # Calculate the T90
        # -- Get the cumulative sum of all productive pores
        df["cum_productive"] = df["productive"].cumsum()
        # -- Find the timepoint (h) at which the cumulative sum >= 90% of the absolute sum
        t90_min = df[df["cum_productive"] >= 0.9 * df["productive"].sum()].index[0]
        pore_activity["t90_h"] = round(t90_min / 60, 1)

        # Add to the db update
        db_update["pore_activity"] = pore_activity

    def parse_minknow_json(self, db_update):
        """Parse useful stuff from the MinKNOW .json report to add to CouchDB"""

        logger.info(f"Parsing report JSON of run {self.run_id}")

        glob_report_json = glob.glob(self.run_dir + "/report*.json")

        if len(glob_report_json) == 0:
            error_message = f"Run {self.run_id} is marked as finished, but missing .json report file."
            logger.error(error_message)
            raise AssertionError(error_message)
        elif len(glob_report_json) > 1:
            error_message = f"Run {self.run_id} is marked as finished, but contains conflicting .json report files."
            logger.error(error_message)
            raise AssertionError(error_message)

        dict_json_report = json.load(open(glob_report_json[0], "r"))

        # Initialize return dict
        parsed_data = {}

        # These sections of the .json can be added as they are
        for section in [
            "software_versions",
            "host",
            "protocol_run_info",
            "user_messages",
        ]:
            parsed_data[section] = dict_json_report[section]

        # Only parse the last acquisition section, which contains the actual sequencing data
        seq_metadata = dict_json_report["acquisitions"][-1]
        seq_metadata_trimmed = {}

        # -- Run info subsection
        seq_metadata_trimmed["acquisition_run_info"] = {}

        seq_metadata_trimmed["acquisition_run_info"]["yield_summary"] = seq_metadata[
            "acquisition_run_info"
        ]["yield_summary"]

        # -- Run output subsection
        seq_metadata_trimmed["acquisition_output"] = []
        for section in seq_metadata["acquisition_output"]:
            if section["type"] in ["AllData", "SplitByBarcode"]:
                seq_metadata_trimmed["acquisition_output"].append(section)

        # -- Read length subseqtion
        seq_metadata_trimmed["read_length_histogram"] = seq_metadata[
            "read_length_histogram"
        ]

        # Add the trimmed acquisition section to the parsed data
        parsed_data["acquisitions"] = []
        parsed_data["acquisitions"].append(seq_metadata_trimmed)

        # Add the parsed data to the db update
        db_update.update(parsed_data)


    def update_transfer_log(self):
        """Update transfer log with run id and date."""
        try:
            with open(self.transfer_log, 'a') as f:
                tsv_writer = csv.writer(f, delimiter='\t')
                tsv_writer.writerow([self.run_id, str(datetime.now())])
                return True
        except IOError:
            logger.warn('Could not update the transfer logfile for run {}. '
                        'Please make sure it gets updated.'.format(self.run_id, self.transfer_log))
            return False

    def archive_run(self):
        """Move directory to nosync."""
        logger.info('Archiving run ' + self.run_id)
        dir_to_move = str(pathlib.Path(self.run_dir).parent)
        project_archive = os.path.join(self.archive_dir, self.experiment_id)
        if not os.path.exists(project_archive):
            os.mkdir(project_archive)
        try:
            print(dir_to_move, project_archive)
            shutil.move(dir_to_move, project_archive)
            logger.info('Successfully archived {}'.format(self.run_id))
            if not os.listdir(self.top_dir):
                logger.info("Project folder {} is empty. Removing it.".format(self.top_dir))
                os.rmdir(self.top_dir)
            else:
                logger.info("Some data is still left in {}. Keeping it.".format(self.top_dir))  # Might be another run for the same project
            return True
        except shutil.Error as e:
            logger.warn('The following error occurred when archiving {}:\n'
                        '{}'.format(self.run_dir, e))
            return False
