import json
import logging
import os
import shutil
from datetime import datetime

from taca.utils import misc
from taca.utils.filesystem import chdir
from taca.utils.statusdb import ElementRunsConnection

logger = logging.getLogger(__name__)


class Run:
    """Defines an Element run"""

    def __init__(self, run_dir, configuration):
        if not os.path.exists(run_dir):
            raise RuntimeError(f"Could not locate run directory {run_dir}")
        self.run_parameters_parsed = False

        self.run_dir = os.path.abspath(run_dir)
        self.CONFIG = configuration

        self.demux_dir = os.path.join(self.run_dir, "Demultiplexing")
        self.final_sequencing_file = os.path.join(self.run_dir, "RunUploaded.json")
        self.demux_stats_file = os.path.join(
            self.demux_dir,
            "RunStats.json",  # Assumes demux is finished when this file is created
        )
        self.run_manifest_file = os.path.join(self.run_dir, "RunManifest.csv")
        self.run_manifest_zip_file = os.path.join(
            self.CONFIG.get("Aviti").get("manifest_zip_location"),
            self.flowcell_id + ".tar.gz",
        )  # TODO: change and add to taca.yaml
        # TODO, need to be real careful when using the flowcell_id as it is manually entered and can mean three different things

        # Instrument generated files
        self.run_parameters_file = os.path.join(self.run_dir, "RunParameters.json")
        self.run_stats_file = os.path.join(self.run_dir, "RunStats.json")
        self.run_manifest_file_from_instrument = os.path.join(
            self.run_dir, "RunManifest.json"
        )
        self.run_uploaded_file = os.path.join(self.run_dir, "RunUploaded.json")

        self.db = ElementRunsConnection(self.CONFIG["statusdb"], dbname="element_runs")

        # Fields to be set by TACA
        self.status = None

        # Fields that will be set when parsing run parameters
        self.run_name = None
        self.run_id = None
        self.side = None
        self.side_letter = None
        self.run_type = None
        self.flowcell_id = None
        self.instrument_name = None
        self.date = None
        self.operator_name = None

    def __str__(self) -> str:
        if self.run_parameters_parsed:
            return f"ElementRun({self.NGI_run_id})"
        else:
            return f"ElementRun({self.run_dir})"

    @property
    def NGI_run_id(self):
        if self.run_parameters_parsed:
            return f"{self.date}_{self.instrument_name}_{self.side_letter}{self.flowcell_id}"
        else:
            raise RuntimeError(f"Run parameters not parsed for run {self.run_dir}")

    def parse_run_parameters(self) -> None:
        """Parse run-information from the RunParameters.json file"""
        try:
            with open(self.run_parameters_file) as json_file:
                run_parameters = json.load(json_file)
        except FileNotFoundError:
            logger.warning(
                f"Run parameters file not found for {self}, might not be ready yet"
            )
            raise

        # Manually entered, but should be side and flowcell id
        self.run_name = run_parameters.get("RunName")

        self.run_id = run_parameters.get(
            "runID"
        )  # Unique hash that we don't really use
        self.side = run_parameters.get("Side")  # SideA or SideB
        self.side_letter = self.side[-1]  # A or B
        self.run_type = run_parameters.get(
            "RunType"
        )  # Sequencing, wash or prime I believe?
        self.flowcell_id = run_parameters.get("FlowcellID")
        self.instrument_name = run_parameters.get("InstrumentName")
        self.date = run_parameters.get("Date")
        self.operator_name = run_parameters.get("OperatorName")
        self.run_parameters_parsed = True

    def to_doc_obj(self):
        # TODO, are we sure what we should do when the RunParameters.json file is missing?

        # Read in all instrument generated files
        instrument_generated_files = {}
        for file in [
            self.run_parameters_file,
            self.run_stats_file,
            self.run_manifest_file_from_instrument,
            self.run_uploaded_file,
        ]:
            if os.path.exists(file):
                with open(file) as json_file:
                    instrument_generated_files[os.path.basename(file)] = json.load(
                        json_file
                    )
            else:
                instrument_generated_files[os.path.basename(file)] = None

        doc_obj = {
            "run_path": self.run_dir,
            "run_status": self.status,
            "NGI_run_id": self.NGI_run_id,
            "instrument_generated_files": instrument_generated_files,
        }

        return doc_obj

    def check_sequencing_status(self):
        if os.path.exists(self.final_sequencing_file):
            with open(self.final_sequencing_file) as json_file:
                sequencing_outcome = json.load(json_file).get("outcome")
            if sequencing_outcome != "OutcomeCompleted":
                return False
            else:
                return True
        else:
            return False

    def get_demultiplexing_status(self):
        if not os.path.exists(self.demux_dir):
            return "not started"
        elif os.path.exists(self.demux_dir) and not os.path.isfile(
            self.demux_stats_file
        ):
            return "ongoing"
        elif os.path.exists(self.demux_dir) and os.path.isfile(self.demux_stats_file):
            return "finished"
        else:
            return "unknown"

    def status_changed(self):
        if not self.run_parameters_parsed:
            raise RuntimeError(
                f"Run parameters not parsed for run {self.run_dir}, cannot check status"
            )
        db_run_status = self.db.check_db_run_status(self.NGI_run_id)
        return db_run_status != self.status

    def update_statusdb(self):
        doc_obj = self.to_doc_obj()
        self.db.upload_to_statusdb(doc_obj)

    def manifest_exists(self):
        return os.path.isfile(self.run_manifest_zip_file)

    def copy_manifests(self):
        shutil.copy(self.run_manifest_zip_file, self.run_dir)
        # TODO: unzip

    def generate_demux_command(self, run_manifest, demux_dir):
        command = [
            self.CONFIG.get(self.software)[
                "bin"
            ],  # TODO add path to bases2fastq executable to config
            self.run_dir,
            demux_dir,
            "-p 12",
        ]
        return command

    def start_demux(self, run_manifest, demux_dir):
        with chdir(self.run_dir):
            cmd = self.generate_demux_command(run_manifest, demux_dir)
            misc.call_external_command_detached(
                cmd, with_log_files=True, prefix="demux_"
            )
            logger.info(
                "Bases2Fastq conversion and demultiplexing "
                f"started for run {os.path.basename(self.run_dir)} on {datetime.now()}"
            )

    def is_transferred(self, transfer_file):
        # TODO: return true if run in transfer log, else false
        pass

    def transfer_ongoing(self):
        # TODO: return true if hidden transfer file marker exists, else false
        pass

    def sync_metadata(self):
        # TODO: copy metadata from demuxed run to ngi-nas-ns
        pass

    def make_transfer_indicator(self):
        # TODO: touch a hidden file in the run directory
        pass

    def transfer(self):
        # TODO: rsync run to analysis cluster
        pass

    def remove_transfer_indicator(self):
        # TODO: remove hidden file in run directory
        pass

    def update_transfer_log(self, transfer_file):
        # TODO: update the transfer log
        pass

    def archive(self):
        # TODO: move run dir to nosync
        pass

    def aggregate_demux_results(self):
        # TODO: aggregate demux results
        pass
