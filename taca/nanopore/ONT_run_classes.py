import os
import logging
import csv
import shutil
import glob
import re
import json
import pandas as pd
import subprocess
import os

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

        # Get run type
        self.run_type = "qc_run" if self.sample_name[0:3] == "QC_" else "user_run"

        # Get attributes from config
        self.minknow_reports_dir = CONFIG["nanopore_analysis"]["minknow_reports_dir"]
        self.analysis_server = CONFIG["nanopore_analysis"]["analysis_server"]
        self.rsync_options = CONFIG["nanopore_analysis"]["rsync_options"]
        self.db = NanoporeRunsConnection(CONFIG["statusdb"], dbname="nanopore_runs")

        self.transfer_details = CONFIG["nanopore_analysis"]["run_types"][self.run_type][
            "instruments"
        ][self.instrument]

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

    def assert_contents(self):
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

    def is_transferred(self) -> bool:
        """Return True if run ID in transfer.tsv, else False."""
        with open(self.transfer_details["transfer_log"], "r") as f:
            return self.run_name in f.read()

    # DB update

    def touch_db_entry(self):
        """Check run vs statusdb. Create entry if there is none."""

        if not self.db.check_run_exists(self):
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

            self.db.create_ongoing_run(self, run_path_file, pore_count_history_file)
            logger.info(
                f"{self.run_name}: Successfully created database entry for ongoing run."
            )
        else:
            logger.info(f"{self.run_name}: Database entry already exists, skipping.")

    def update_db_entry(self, force_update=False):
        """Check run vs statusdb. Create or update run entry."""

        # If no run document exists in the database, ceate an ongoing run document
        self.touch_db_entry()

        # If the run document is marked as "ongoing" or database is being manually updated
        if self.db.check_run_status(self) == "ongoing" or force_update == True:
            logger.info(
                f"{self.run_name}: Run exists in the database with run status: {self.db.check_run_status(self)}."
            )

            logger.info(f"{self.run_name}: Updating...")

            # Instantiate json (dict) to update the db with
            db_update = {}

            # Parse report_*.json
            self.parse_minknow_json(db_update)

            # Parse pore_activity_*.csv
            self.parse_pore_activity(db_update)

            # Update the DB entry
            self.db.finish_ongoing_run(self, db_update)

        # If the run document is marked as "finished"
        elif self.db.check_run_status(self) == "finished":
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

    def copy_metadata(self):
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
        dst = self.transfer_details["metadata_dir"]

        os.system(
            f"rsync -rv --exclude={{{','.join(exclude_patterns_quoted)}}} {src} {dst}"
        )

    def copy_html_report(self):

        logger.info(f"{self.run_name}: Transferring .html report to ngi-internal...")

        # Transfer the MinKNOW .html report file to ngi-internal, renaming it to the full run ID. Requires password-free SSH access.
        report_src_path = self.get_file("/report*.html")
        report_dest_path = os.path.join(
            self.minknow_reports_dir,
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

    # Transfer run

    def transfer_run(self):
        """Transfer dir to destination specified in config file via rsync"""
        destination = self.transfer_details["destination"]

        for k, v in self.rsync_options.items():
            if v == "None":
                self.rsync_options[k] = None
        logger.info(
            f"{self.run_name}: Transferring to {self.analysis_server['host'] if self.analysis_server else destination}..."
        )
        if self.analysis_server:
            transfer_object = RsyncAgent(
                self.run_abspath,
                dest_path=destination,
                remote_host=self.analysis_server["host"],
                remote_user=self.analysis_server["user"],
                validate=False,
                opts=self.rsync_options,
            )
        else:
            transfer_object = RsyncAgent(
                self.run_abspath,
                dest_path=destination,
                validate=False,
                opts=self.rsync_options,
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
            with open(self.transfer_details["transfer_log"], "a") as f:
                tsv_writer = csv.writer(f, delimiter="\t")
                tsv_writer.writerow([self.run_name, str(datetime.now())])
        except IOError:
            msg = f"{self.run_name}: Could not update the transfer logfile {self.transfer_details['transfer_log']}"
            logger.error(msg)
            raise IOError(msg)

    # Archive run

    def archive_run(self):
        """Move directory to nosync."""
        logger.info(f"{self.run_name}: Archiving run...")

        src = self.run_abspath
        dst = os.path.join(self.run_abspath, os.pardir, "nosync")

        shutil.move(src, dst)
        logger.info(f"{self.run_name}: Archiving run successful.")


class ONT_user_run(ONT_run):
    """ONT user run, has class methods and attributes specific to user runs."""

    def __init__(self, run_abspath: str):
        super(ONT_user_run, self).__init__(run_abspath)


class ONT_qc_run(ONT_run):
    """ONT QC run, has class methods and attributes specific to QC runs"""

    def __init__(self, run_abspath: str):
        super(ONT_qc_run, self).__init__(run_abspath)

        # Get Anglerfish attributes from run
        self.anglerfish_done_abspath = f"{self.run_abspath}/.anglerfish_done"
        self.anglerfish_ongoing_abspath = f"{self.run_abspath}/.anglerfish_ongoing"

        # Get Anglerfish attributes from config
        self.anglerfish_config = CONFIG["nanopore_analysis"]["run_types"][
            self.run_type
        ]["anglerfish"]

        self.anglerfish_samplesheets_dir = self.anglerfish_config[
            "anglerfish_samplesheets_dir"
        ]
        self.anglerfish_samplesheets_dir = self.anglerfish_config[
            "anglerfish_samplesheets_dir"
        ]
        self.anglerfish_env_name = self.anglerfish_config["anglerfish_env_name"]
        # Depending on the Conda installation, Conda may need to be initialized when used in a subprocess
        if "conda_init_path" in self.anglerfish_config:
            self.conda_init_path = self.anglerfish_config["conda_init_path"]
        else:
            self.conda_init_path = None

    # QC methods

    def get_anglerfish_exit_code(self) -> int or None:
        """Check whether Anglerfish has finished.

        Return exit code or None.
        """
        if os.path.exists(self.anglerfish_done_abspath):
            return int(open(self.anglerfish_done_abspath, "r").read())
        else:
            return None

    def get_anglerfish_pid(self) -> str or None:
        """Check whether Anglerfish is ongoing.

        Return process ID or None."""
        if os.path.exists(self.anglerfish_ongoing_abspath):
            return str(open(self.anglerfish_ongoing_abspath, "r").read())
        else:
            return None

    def fetch_anglerfish_samplesheet(self) -> str:
        """Fetch Anglerfish samplesheet belonging to the run from where it was
        dumped by LIMS and put it in the run directory.

        a) If file(s) is available, copy to run folder. On success, return True and
        add new samplesheet abspath as run object attribute.

        b) If the file is not yet available, return False.
        """

        # Following line assumes run was started same year as samplesheet was generated
        current_year = self.date[0:4]
        expected_file_pattern = f"anglerfish_samplesheet_{self.experiment_name}_*.csv"

        # Finalize query pattern
        pattern_abspath = os.path.join(
            self.anglerfish_samplesheets_dir, current_year, expected_file_pattern
        )

        glob_results = glob.glob(pattern_abspath)

        if len(glob_results) == 0:
            return False

        else:
            # Sort by ascending date
            glob_results.sort()

            # Grab abspath of latest samplesheet
            src = glob_results[-1]
            dst = os.path.join(self.run_abspath, os.path.basename(src))

            # Copy into run directory
            if os.system(f"rsync -v {src} {dst}") == 0:
                self.anglerfish_samplesheet = dst
                return True
            else:
                raise RsyncError(
                    f"{self.run_name}: Error occured when copying anglerfish samplesheet to run dir."
                )

    def run_anglerfish(self):
        """Run Anglerfish as subprocess within it's own Conda environment.
        Dump files to indicate ongoing and finished processes.
        """

        ss_basename = os.path.basename(self.anglerfish_samplesheet)
        pid_yymmdd_hhmm = ss_basename.split(".")[0].split("_")[-3:]
        anglerfish_run_name = f"anglerfish_run_using_ss_{'_'.join(pid_yymmdd_hhmm)}_on"

        n_threads = 2  # This could possibly be changed

        anglerfish_command = (
            "anglerfish"
            + f" --samplesheet {self.anglerfish_samplesheet}"
            + f" --out_fastq {self.run_abspath}"
            + f" --run_name {anglerfish_run_name}"
            + f" --threads {n_threads}"
            + f" --skip_demux"
        )

        full_command = (
            # Some systems need Conda to be initialized in the subshell
            (f"source {self.conda_init_path} && " if self.conda_init_path else "")
            # Activate environment
            + f"conda activate {self.anglerfish_env_name}"
            # On success: Run Anglerfish
            + f" && {anglerfish_command}"
            # Dump Anglerfish exit status in file
            + f" ; echo $? > {self.anglerfish_done_abspath}"
            # Remove file indicating ongoing run
            + f" ; rm {self.anglerfish_ongoing_abspath}"
        )

        try:
            # Start Anglerfish
            with subprocess.Popen(
                full_command,
                shell=True,
                stdout=subprocess.PIPE,
                encoding="utf-8",
                cwd=self.run_abspath,
            ) as process:
                # Create file indicating ongoing run, containing process id
                os.system(f"echo '{process.pid}' > {self.anglerfish_ongoing_abspath}")
                logger.info(
                    f"{self.run_name}: Anglerfish started with process ID {process.pid}."
                )

        except subprocess.CalledProcessError:
            logger.warn(f"{self.run_name}: An error occured when running Anglerfish.")
