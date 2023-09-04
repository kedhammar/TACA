import os
import subprocess
import shutil
import glob
import logging
import pathlib

from taca.nanopore.ONT_run import ONT_run
from taca.utils.config import CONFIG
from taca.utils.minion_barcodes import BARCODES

logger = logging.getLogger(__name__)


class MinIONqc(ONT_run):
    """Minion QC run"""
    def __init__(self, run_dir, anglerfish_sample_sheet):
        super(MinIONqc, self).__init__(run_dir)
        self.transfer_details = (
            CONFIG.get("nanopore_analysis").get("minion_qc_run").get("transfer")
        )
        self.transfer_log = self.transfer_details.get("transfer_file")
        self.archive_dir = self.transfer_details.get("finished_dir")

        self.anglerfish_dir = os.path.join(self.run_abspath, "anglerfish_output")
        self.anglerfish_exit_status_file = os.path.join(
            self.run_abspath, ".exitcode_for_anglerfish"
        )
        self.year_processed = self.run_name[0:4]
        self.flowcell_id = self.run_name.split("_")[3]
        self.anglerfish_sample_sheet = anglerfish_sample_sheet

    def get_anglerfish_samplesheet(self):
        """Copy Anglerfish sample sheet from LIMS."""
        lims_samplesheet_dir = os.path.join(
            CONFIG.get("nanopore_analysis")
            .get("minion_qc_run")
            .get("samplesheets_dir"),
            self.year_processed,
        )
        found_samplesheets = glob.glob(
            lims_samplesheet_dir + "/*" + self.experiment_name + "*"
        )
        if not found_samplesheets:
            logger.warn(
                "No Anglerfish sample sheet from LIMS found for run {}.".format(
                    self.run_name
                )
            )
            return None
        elif len(found_samplesheets) > 1:
            logger.warn(
                "Found more than one Anglerfish sample sheet from LIMS for run {}.".format(
                    self.run_name
                )
            )
            return None
        else:
            sample_sheet_copy = os.path.join(
                self.run_abspath, os.path.basename(found_samplesheets[0])
            )
            shutil.copyfile(
                found_samplesheets[0], sample_sheet_copy
            )  # This will overwrite any existing file with the same name. If you want to use a manually edited sample sheet, name it something else and run TACA manually
            return sample_sheet_copy

    def check_exit_status(self, status_file):
        """Read pipeline exit status file and return True if 0, False if anything else"""
        with open(status_file, 'r') as f:
            exit_status = f.readline().strip()
        return exit_status == '0'

    def start_anglerfish(self):
        """Start Anglerfish."""
        os.makedirs(self.anglerfish_dir)
        anglerfish_command = ('anglerfish'
                            + ' --samplesheet ' + self.anglerfish_sample_sheet
                            + ' --out_fastq ' + self.anglerfish_dir
                            + ' --threads 2'
                            + ' --skip_demux; echo $? > .exitcode_for_anglerfish')
        try:
            p_handle = subprocess.Popen(
                anglerfish_command,
                stdout=subprocess.PIPE,
                shell=True,
                cwd=self.run_abspath,
            )
            logger.info(
                "Started Anglerfish for run {} using: {}".format(
                    self.run_abspath, anglerfish_command
                )
            )
        except subprocess.CalledProcessError:
            logger.warn(
                "An error occurred while starting the Anglerfish for run {}. "
                "Please check the logfile for info.".format(self.run_abspath)
            )

    def copy_results_for_lims(self):
        """Find results and copy to lims directory."""
        lims_result_path = os.path.join(
            CONFIG.get("nanopore_analysis")
            .get("minion_qc_run")
            .get("lims_results_dir"),
            self.run_name,
        )
        lims_result_file = os.path.join(
            lims_result_path, "anglerfish_stats_" + self.experiment_name + ".txt"
        )
        anglerfish_results = self._find_anglerfish_results()
        if not os.path.isdir(lims_result_path):
            os.mkdir(lims_result_path)
        try:
            shutil.copyfile(anglerfish_results, lims_result_file)
            return True
        except TypeError as e:
            logger.warn(
                "An error occurred while copying the Anglerfish results for {} to lims: {}".format(
                    self.run_name, e
                )
            )
            return False

    def _find_anglerfish_results(self):
        """Return location of Anglerfish results."""
        results_file = ''
        for sub_dir in os.listdir(self.anglerfish_dir):
            if 'anglerfish_stats.txt' in os.listdir(os.path.join(self.anglerfish_dir, sub_dir)):
                results_file = os.path.join(self.anglerfish_dir, sub_dir, 'anglerfish_stats.txt')
                return results_file
        if not results_file:
            logger.warn('Could not find any Anglerfish results in {}'.format(self.anglerfish_dir))

class MinIONdelivery(ONT_run):
    """Minion delivery run"""
    def __init__(self, run_dir):
        super(MinIONdelivery, self).__init__(run_dir)
        self.transfer_details = CONFIG.get('nanopore_analysis').get('minion_delivery_run').get('transfer')
        self.transfer_log = self.transfer_details.get('transfer_file')
        self.archive_dir = self.transfer_details.get('finished_dir')
    
    def dump_path(self):
        """Dump path to run to a file that can be
        used when uploading stats to statusdb from preproc."""
        new_file = os.path.join(self.run_abspath, "run_path.txt")
        proj, sample, run = self.run_abspath.split("/")[-3:]
        path_to_write = os.path.join(proj, sample, run)
        with open(new_file, 'w') as f:
            f.write(path_to_write)
    
    def write_finished_indicator(self):
        """Write a hidden file to indicate 
        when the finial rsync is finished."""
        new_file = os.path.join(self.run_abspath, ".sync_finished")
        pathlib.Path(new_file).touch()
        return new_file