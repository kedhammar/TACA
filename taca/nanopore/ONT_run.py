import os
import logging
import csv
import shutil
import glob
import re
import json
import pandas as pd

from taca.utils.statusdb import NanoporeRunsConnection
from datetime import datetime
from taca.utils.config import CONFIG
from taca.utils.transfer import RsyncAgent, RsyncError

logger = logging.getLogger(__name__)

ONT_RUN_PATTERN = re.compile(
    "^(\d{8})_(\d{4})_([0-9a-zA-Z]+)_([0-9a-zA-Z]+)_([0-9a-zA-Z]+)$"
)


class ONT_run(object):
    """General Nanopore run.

    Expects instantiation from absolute path of run directory on preprocessing server.
    """

    def __init__(self, run_abspath: str):

        # Get paths and names of MinKNOW experiment, sample and run

        self.run_name = os.path.basename(run_abspath)
        self.run_abspath = run_abspath

        assert re.match(
            ONT_RUN_PATTERN, self.run_name
        ), f"Run {self.run_name} doesn't look like a run dir"

        # Parse MinKNOW sample and experiment name
        with open(self.get_file("/run_path.txt"), "r") as stream:
            self.experiment_name, self.sample_name, _ = stream.read().split("/")

        # Get info from run name
        (
            self.date,
            self.time,
            self.position,
            self.flowcell_id,
            self.run_hash,
        ) = self.run_name.split("_")

        # Get instrument
        self.instrument = "promethion" if len(self.position) == 2 else "minion"

        # Get instument-specific transfer details
        self.transfer_details = (
            CONFIG.get("nanopore_analysis").get("ont_transfer").get(self.instrument)
        )
        self.transfer_log = self.transfer_details.get("transfer_file")
        self.archive_dir = self.transfer_details.get("finished_dir")
        self.metadata_dir = self.transfer_details.get("metadata_dir")

        # Get QC bool
        self.qc = True if self.sample_name[0:3] == "QC_" else False

    # Looking for files within the run dir

    def has_file(self, content_pattern: str) -> bool:
        """Checks within run dir for pattern, e.g. '/report*.json', returns bool."""
        query_path = self.run_abspath + content_pattern
        query_glob = glob.glob(query_path)

        if len(query_glob) > 0:
            return True
        else:
            return False

    def get_file(self, content_pattern) -> str:
        """Checks within run dir for pattern, e.g. '/report*.json', returns file abspath as string."""
        query_path = self.run_abspath + content_pattern
        query_glob = glob.glob(query_path)

        if len(query_glob) == 1:
            return query_glob[0]
        elif len(query_glob) == 0:
            raise AssertionError(f"Could not find {query_path}")
        else:
            raise AssertionError(f"Found multiple instances of {query_path}")

    # Evaluating run status

    def is_synced(self) -> bool:
        return self.has_file("/.sync_finished")

    def is_transferred(self) -> bool:
        """Return True if run ID in transfer.tsv, else False."""
        with open(self.transfer_log, "r") as f:
            return self.run_name in f.read()

    def assert_contents(self) -> bool:
        """Checklist function to assure run has all files necessary to proceed with processing"""

        # Completion indicators
        assert self.has_file("/.sync_finished")
        assert self.has_file("/final_summary*.txt")

        # NGI files from instrument
        assert self.has_file("/pore_count_history.csv")
        assert self.has_file("/run_path.txt")

        # MinKNOW reports
        assert self.has_file("/report_*.json")
        assert self.has_file("/report_*.html")

        # MinKNOW auxillary files
        assert self.has_file("/final_summary*.txt")
        assert self.has_file("/pore_activity*.csv")

    # DB update

    def touch_db_entry(self):
        """Check run vs statusdb. Create entry if there is none."""

        sesh = NanoporeRunsConnection(CONFIG["statusdb"], dbname="nanopore_runs")

        if not sesh.check_run_exists(self):
            logger.info(
                f"{self.run_name}: Run does not exist in the database, creating entry for ongoing run."
            )

            run_path_file = os.path.join(self.run_abspath, "run_path.txt")
            assert os.path.isfile(run_path_file), f"Couldn't find {run_path_file}"

            pore_count_history_file = os.path.join(
                self.run_abspath, "pore_count_history.csv"
            )
            assert os.path.isfile(
                pore_count_history_file
            ), f"Couldn't find {pore_count_history_file}"

            sesh.create_ongoing_run(self, run_path_file, pore_count_history_file)
            logger.info(
                f"{self.run_name}: Successfully created database entry for ongoing run."
            )
        else:
            logger.info(f"{self.run_name}: Database entry already exists, skipping.")

    def update_db_entry(self, force_update=False):
        """Check run vs statusdb. Create or update run entry."""

        sesh = NanoporeRunsConnection(CONFIG["statusdb"], dbname="nanopore_runs")

        # If no run document exists in the database, ceate an ongoing run document
        self.touch_db_entry()

        # If the run document is marked as "ongoing" or database is being manually updated
        if sesh.check_run_status(self) == "ongoing" or force_update == True:
            logger.info(
                f"{self.run_name}: Run exists in the database with run status: {sesh.check_run_status(self)}."
            )

            logger.info(f"{self.run_name}: Updating...")

            # Instantiate json (dict) to update the db with
            db_update = {}

            # Parse report_*.json
            self.parse_minknow_json(db_update)

            # Parse pore_activity_*.csv
            self.parse_pore_activity(db_update)

            # Update the DB entry
            sesh.finish_ongoing_run(self, db_update)

        # If the run document is marked as "finished"
        elif sesh.check_run_status(self) == "finished":
            logger.info(
                f"Run {self.run_name} exists in the database as an finished run, do nothing."
            )

    def parse_pore_activity(self, db_update):

        logger.info(f"{self.run_name}: Parsing pore activity...")

        pore_activity = {}

        # Use pandas to pivot the data into a more manipulable dataframe
        df = pd.read_csv(self.get_file("/pore_activity_*.csv"))
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

        logger.info(f"{self.run_name}:Parsing report JSON...")

        dict_json_report = json.load(open(self.get_file("/report*.json"), "r"))

        # Initialize return dict
        parsed_data = {}

        # These sections of the .json can be added as they are
        for section in [
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

    # Transferring metadata

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

        src = self.run_abspath
        dst = self.metadata_dir

        os.system(
            f"rsync -rv --exclude={{{','.join(exclude_patterns_quoted)}}} {src} {dst}"
        )

    def transfer_html_report(self):

        logger.info(f"{self.run_name}: Transferring .html report to ngi-internal...")

        # Transfer the MinKNOW .html report file to ngi-internal, renaming it to the full run ID. Requires password-free SSH access.
        report_src_path = self.get_file("/report*.html")
        report_dest_path = os.path.join(
            CONFIG["nanopore_analysis"]["ont_transfer"]["minknow_reports_dir"],
            f"report_{self.run_name}.html",
        )
        transfer_object = RsyncAgent(
            src_path=report_src_path,
            dest_path=report_dest_path,
            validate=False,
        )
        try:
            transfer_object.transfer()
            logger.info(
                f"Successfully transferred the MinKNOW report of run {self.run_name}"
            )
        except RsyncError:
            msg = f"{self.run_name}: An error occurred while attempting to transfer the report {report_src_path} to {report_dest_path}."
            logger.error(msg)
            raise RsyncError(msg)

    # Transfer data

    def transfer_run(self) -> bool:
        """Transfer dir to destination specified in config file via rsync"""
        destination = self.transfer_details.get("destination")
        rsync_opts = self.transfer_details.get("rsync_options")
        for k, v in rsync_opts.items():
            if v == "None":
                rsync_opts[k] = None
        connection_details = self.transfer_details.get("analysis_server", None)
        logger.info(
            f"{self.run_name}: Transferring to {connection_details['host'] if connection_details else destination}..."
        )
        if connection_details:
            transfer_object = RsyncAgent(
                self.run_abspath,
                dest_path=destination,
                remote_host=connection_details.get("host"),
                remote_user=connection_details.get("user"),
                validate=False,
                opts=rsync_opts,
            )
        else:
            transfer_object = RsyncAgent(
                self.run_abspath, dest_path=destination, validate=False, opts=rsync_opts
            )
        try:
            transfer_object.transfer()
        except RsyncError:
            msg = f"{self.run_name}: An error occurred while transferring to the analysis server."
            logger.error(msg)
            raise RsyncError(msg)

    def update_transfer_log(self):
        """Update transfer log with run id and date."""
        try:
            with open(self.transfer_log, "a") as f:
                tsv_writer = csv.writer(f, delimiter="\t")
                tsv_writer.writerow([self.run_name, str(datetime.now())])
        except IOError:
            msg = f"{self.run_name}: Could not update the transfer logfile {self.transfer_log}"
            logger.error(msg)
            raise IOError(msg)

    def archive_run(self):
        """Move directory to nosync."""
        logger.info(f"{self.run_name}: Archiving run...")

        src = self.run_abspath
        dst = os.path.join(self.run_abspath, os.pardir, "nosync")

        shutil.move(src, dst)
        logger.info(f"{self.run_name}: Archiving run successful.")
