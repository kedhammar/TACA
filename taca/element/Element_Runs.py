import json
import logging
import os
import re
import csv
import zipfile
import subprocess
import shutil
from datetime import datetime
from pathlib import Path
from glob import glob

import pandas as pd

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
        self.demux_stats_file = "RunStats.json"  # Assumes demux is finished when this file is created
        self.transfer_file = (
            self.CONFIG.get("Element").get(self.sequencer_type).get("transfer_log")
        )  # TODO: change and add to taca.yaml
        self.rsync_exit_file = os.path.join(self.run_dir, '.rsync_exit_status')

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
        self.lims_step_id = None
        self.lims_full_manifest = None
        self.lims_start_manifest = None
        self.lims_demux_manifests = None

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
        # TODO: are we sure what we should do when the RunParameters.json file is missing?

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
        demux_dirs = glob.glob(
            os.path.join(self.run_dir, "Delmultiplexing*")
            )
        finished_count = 0
        for demux_dir in demux_dirs:
            if os.path.exists(self.demux_dir) and not os.path.isfile(
                os.path.join(demux_dir, self.demux_stats_file)
                ):
                return "ongoing"
            elif os.path.exists(self.demux_dir) and os.path.isfile(
                os.path.join(demux_dir, self.demux_stats_file)
                ):
                finished_count += 1  # TODO: check exit status of demux in exit status file
        if finished_count == len(demux_dirs):
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

    def get_lims_step_id(self) -> str | None:
        """If the run was started using a LIMS-generated manifest,
        the ID of the LIMS step can be extracted from it.
        """

        # TODO test me

        assert self.manifest_exists(), "Run manifest not found"
        with open(self.run_manifest_file_from_instrument) as csv_file:
            manifest_lines = csv_file.readlines()
        for line in manifest_lines:
            if "lims_step_id" in line:
                lims_step_id = line.split(",")[1]
                return lims_step_id
        return None

    def copy_manifests(self) -> bool:
        """Fetch the LIMS-generated run manifests from ngi-nas-ns and unzip them into a run subdir."""

        # TODO test me

        # Specify dir in which LIMS drop the manifest zip files
        dir_to_search = os.path.join(
            self.CONFIG.get("Aviti").get(
                "manifest_zip_location"
            ),  # TODO: change and add to taca.yaml
            datetime.now().year,
        )

        # Use LIMS step ID if available, else flowcell ID, to make a query pattern
        if self.lims_step_id:
            logging.info(
                f"Using LIMS step ID '{self.lims_step_id}' to find LIMS run manifests."
            )
            glob_pattern = f"{dir_to_search}/*{self.lims_step_id}*.zip"
        else:
            logging.warning(
                "LIMS step ID not available, using flowcell ID to find LIMS run manifests."
            )
            glob_pattern = f"{dir_to_search}/*{self.flowcell_id}*.zip"

        # Find paths matching the pattern
        glob_results = glob(glob_pattern)
        if len(glob_results) == 0:
            logger.warning(
                f"No manifest found for run '{self.run_dir}' with pattern '{glob_pattern}'."
            )
            return False  # TODO determine whether to raise an error here instead
        elif len(glob_results) > 1:
            logger.warning(
                f"Multiple manifests found for run '{self.run_dir}' with pattern '{glob_pattern}', using latest one."
            )
            glob_results.sort()
            zip_src_path = glob_results[-1]
        else:
            zip_src_path = glob_results[0]

        # Make a run subdir named after the zip file and extract manifests there
        zip_name = os.path.basename(zip_src_path)
        zip_dst_path = os.path.join(self.run_dir, zip_name)
        os.mkdir(zip_dst_path)

        with zipfile.ZipFile(zip_src_path, "r") as zip_ref:
            zip_ref.extractall(zip_dst_path)

        # Set the paths of the different manifests as attributes
        manifests = os.listdir(zip_dst_path)
        self.lims_full_manifest = [
            m for m in manifests if re.match(r".*_untrimmed\.csv$", m)
        ][0]
        self.lims_start_manifest = [
            m for m in manifests if re.match(r".*_trimmed\.csv$", m)
        ][0]
        self.lims_demux_manifests = [
            m for m in manifests if re.match(r".*_\d+\.csv$", m)
        ]

        return True

    def make_demux_manifests(
        self, manifest_to_split: os.PathLike, outdir: os.PathLike | None = None
    ) -> list[os.PathLike]:
        """Derive composite demultiplexing manifests (grouped by index duplicity and lengths)
        from a single information-rich manifest.
        """

        # TODO test me

        # Read specified manifest
        with open(manifest_to_split) as f:
            manifest_contents = f.read()

        # Get '[SAMPLES]' section
        split_contents = "[SAMPLES]".split(manifest_contents)
        assert (
            len(split_contents) == 2
        ), f"Could not split sample rows out of manifest {manifest_contents}"
        sample_section = split_contents[1].split("\n")

        # Split into header and rows
        header = sample_section[0]
        sample_rows = sample_section[1:]

        # Convert to list of dicts
        sample_dicts = []
        for row in sample_rows:
            row_dict = dict(zip(header.split(","), row.split(",")))
            sample_dicts.append(row_dict)

        # Convert to dataframe
        df = pd.DataFrame.from_dict(sample_dicts)

        # Separate samples from controls
        df_samples = df[df["Project"] != "Control"].copy()
        df_controls = df[df["Project"] == "Control"].copy()

        # Apply default dir path for output
        if outdir is None:
            outdir = self.run_dir

        ## Build composite manifests

        manifest_root_name = f"{self.NGI_run_id}_demux"

        # Get idx lengths for calculations
        df_samples.loc[:, "len_idx1"] = df["Index1"].apply(len)
        df_samples.loc[:, "len_idx2"] = df["Index2"].apply(len)

        # Break down by index lengths and lane, creating composite manifests
        manifests = []
        n = 0
        for (len_idx1, len_idx2, lane), group in df_samples.groupby(
            ["len_idx1", "len_idx2", "Lane"]
        ):
            file_name = f"{manifest_root_name}_{n}.csv"
            runValues_section = "\n".join(
                [
                    "[RUNVALUES]",
                    "KeyName, Value",
                    f'manifest_file, "{file_name}"',
                    f"manifest_group, {n+1}/{len(df.groupby(['len_idx1', 'len_idx2', 'Lane']))}",
                    f"grouped_by, len_idx1:{len_idx1} len_idx2:{len_idx2} lane:{lane}",
                ]
            )

            settings_section = "\n".join(
                [
                    "[SETTINGS]",
                    "SettingName, Value",
                ]
            )

            # Add PhiX stratified by index length
            if group["phix_loaded"].any():
                # Subset controls by lane
                group_controls = df_controls[df_controls["Lane"] == lane].copy()

                # Trim PhiX indexes to match group
                group_controls.loc[:, "Index1"] = group_controls.loc[:, "Index1"].apply(
                    lambda x: x[:len_idx1]
                )
                group_controls.loc[:, "Index2"] = group_controls.loc[:, "Index2"].apply(
                    lambda x: x[:len_idx2]
                )

                # Add PhiX to group
                group = pd.concat([group, group_controls], axis=0, ignore_index=True)

            samples_section = (
                f"[SAMPLES]\n{group.iloc[:, 0:6].to_csv(index=None, header=True)}"
            )

            manifest_contents = "\n\n".join(
                [runValues_section, settings_section, samples_section]
            )

            file_path = os.path.join(outdir, file_name)
            manifests.append((file_path, manifest_contents))
            n += 1

        for manifest_path, manifest_contents in manifests:
            with open(os.path.join(outdir, manifest_path), "w") as f:
                f.write(manifest_contents)

        manifest_paths = [t[0] for t in manifests]
        return manifest_paths

    def generate_demux_command(self, run_manifest, demux_dir):
        command = (f"{self.CONFIG.get(self.software)["bin"]}"   # TODO: add path to bases2fastq executable to config
            + f" {self.run_dir}"
            + f" {demux_dir}"
            + " -p 8"
            + f" -r {run_manifest}"
            + " --legacy-fastq"  # TODO: except if Smart-seq3
            + f" --force-index-orientation"
            )  # TODO: any other options?
        return command

    def start_demux(self, run_manifest, demux_dir):
        with chdir(self.run_dir):
            cmd = self.generate_demux_command(run_manifest, demux_dir)
            # TODO: handle multiple composite manifests for demux
            try:
                p_handle = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, cwd=self.run_dir)
                logger.info(
                    "Bases2Fastq conversion and demultiplexing "
                    f"started for run {self} on {datetime.now()}"
                )
            except subprocess.CalledProcessError:
                logger.warning("An error occurred while starting demultiplexing for "
                               f"{self} on {datetime.now()}."
                )
        return


    def get_transfer_status(self):
        if not self.in_transfer_log() and not self.transfer_ongoing() and not self.rsync_complete():
            return "not started"
        elif self.transfer_ongoing() and not self.rsync_complete():
            return "ongoing"
        elif self.rsync_complete() and not self.in_transfer_log():
            return "rsync done"
        elif self.in_transfer_log():
            return "unknown"

    def in_transfer_log(self):
        with open(self.transfer_file, "r") as transfer_file:
            for row in transfer_file.read():
                if self.NGI_run_id in row:
                    return True
        return False

    def transfer_ongoing(self):
        return os.path.isfile(os.path.join(self.run_dir, '.rsync_ongoing'))

    def rsync_complete(self):
        return os.path.isfile(self.rsync_exit_file)

    def rsync_successful(self):
        with open(os.path.join(self.run_dir, '.rsync_exit_status')) as rsync_exit_file:
            rsync_exit_status = rsync_exit_file.readlines()
        if rsync_exit_status[0].strip() == 0:
            return True
        else:
            return False

    # Clear all content under a dir
    def clear_dir(dir):
        for filename in os.listdir(dir):
            file_path = os.path.join(dir, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path} Reason {e}")

    # Create symlink for a simple demultiplexing dir
    def symlink_demux_dir(src_dir, dest_dir):
        # Ensure the destination directory exists
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
        # Clear all content under dest_dir
        clear_dir(dest_dir)
        # Loop through all files and directories in the source directory
        for item in os.listdir(src_dir):
            src_path = os.path.join(src_dir, item)
            # Move content of Samples to the parental dir
            if item == "Samples":
                dest_path = dest_dir
            else:
                dest_path = os.path.join(dest_dir, item)
            try:
                # Create symbolic link only if it doesn't already exist
                if not os.path.exists(dest_path):
                    os.symlink(src_path, dest_path)
                    print(f"Linked {src_path} to {dest_path}")
                else:
                    print(f"{dest_path} already exists.")
            except OSError as e:
                print(f"Error linking {src_path} to {dest_path}: {e}")


    # Collect demux info into a list of dictionaries
    # Structure: [{'sub_demux_count':XXX, 'SampleName':XXX, 'Index1':XXX, 'Index2':XXX, 'Lane':XXX, 'Project':XXX, 'Recipe':XXX}]
    def collect_demux_runmanifest(self, demux_results_dirs):
        demux_runmanifest = []
        for demux_dir in demux_results_dirs:
            sub_demux_count = demux_dir.split('_')[1]
            with open(os.path.join(self.run_dir, demux_dir, 'RunManifest.csv'), 'r') as file:
                lines = file.readlines()
            sample_section = False
            headers = []
            # Loop through each line
            for line in lines:
                # Check if we reached the "[SAMPLES]" section
                if '[SAMPLES]' in line:
                    sample_section = True
                    continue
                # Exit the sample section if another section is encountered
                if sample_section and line.startswith('['):
                    break
                # If in the sample section, process the sample lines
                if sample_section:
                    # Clean up the line
                    line = line.strip()
                    # Skip empty lines
                    if not line:
                        continue
                    # Get the headers from the first line
                    if not headers:
                        headers = line.split(',')
                    else:
                        # Parse sample data
                        values = line.split(',')
                        sample_dict = dict(zip(headers, values))
                        sample_dict['sub_demux_count'] = sub_demux_count
                        demux_runmanifest.append(sample_dict)
        sorted_demux_runmanifest = sorted(demux_runmanifest, key=lambda x: (x['Lane'], x['SampleName'], x['sub_demux_count']))
        return sorted_demux_runmanifest


    # Aggregate the output FastQ files of samples from multiple demux
    def aggregate_sample_fastq(self, demux_runmanifest):
        lanes = sorted(list(set(sample['Lane'] for sample in demux_runmanifest)))
        unique_sample_demux = set()
        for lane in lanes:
            sample_count = 1
            for sample in demux_runmanifest:
                lanenr = sample['Lane']
                project = sample['Project']
                sample = sample['SampleName']
                sub_demux_count = sample['sub_demux_count']
                # Skip PhiX
                if lanenr == lane and sample != "PhiX":
                    sample_tuple = (sample, sub_demux_count)
                    if sample_tuple not in unique_sample_demux:
                        project_dest = os.path.join(self.run_dir, self.demux_dir, project)
                        sample_dest = os.path.join(self.run_dir, self.demux_dir, project, sample)
                        if not os.path.exists(project_dest):
                            os.makedirs(project_dest)
                        if not os.path.exists(sample_dest):
                            os.makedirs(sample_dest)
                        fastqfiles = glob.glob(os.path.join(self.run_dir, f"Demultiplexing_{sub_demux_count}", "Samples", project, sample, f"*L00{lane}*.fastq.gz"))
                        for fastqfile in fastqfiles:
                            old_name = os.path.basename(fastqfile)
                            read_label = re.search(rf"L00{lane}_(.*?)_001", old_name).group(1)
                            new_name = "_".join([sample, f"S{sample_count}", f"L00{lane}", read_label, "001.fastq.gz"])
                            os.symlink(fastqfile, os.path.join(sample_dest, new_name))
                        unique_sample_demux.add(sample_tuple)
                        sample_count += 1


    # Symplink the output FastQ files of undet only if a lane does not have multiple demux
    def aggregate_undet_fastq(self, demux_runmanifest):
        lanes = sorted(list(set(sample['Lane'] for sample in demux_runmanifest)))
        for lane in lanes:
            sub_demux = list(set(sample['sub_demux_count'] for sample in demux_runmanifest if sample['Lane']==lane))
            if len(sub_demux) == 1:
                project_dest = os.path.join(self.run_dir, self.demux_dir, "Undetermined")
                if not os.path.exists(project_dest):
                    os.makedirs(project_dest)
                fastqfiles = glob.glob(os.path.join(self.run_dir, f"Demultiplexing_{sub_demux[0]}", "Samples", "Undetermined", "*.fastq.gz"))
                for fastqfile in fastqfiles:
                    base_name = os.path.basename(fastqfile)
                    os.symlink(fastqfile, os.path.join(project_dest, base_name))

    # Write to csv
    def write_to_csv(data, filename):
        # Get the fieldnames from the keys of the first dictionary
        fieldnames = data[0].keys()
        # Open the file and write the CSV
        with open(filename, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            # Write the header (fieldnames)
            writer.writeheader()
            # Write the data (rows)
            writer.writerows(data)


    # Aggregate stats in UnassignedSequences.csv
    def aggregate_stats_unassigned(self, demux_runmanifest):
        aggregated_unassigned_indexes = []
        lanes = sorted(list(set(sample['Lane'] for sample in demux_runmanifest)))
        for lane in lanes:
            sub_demux_index_lens = set()
            for sample in demux_runmanifest:
                if sample['Lane'] == lane:
                    sub_demux_index_lens.add((sample['sub_demux_count'], (len(sample.get("Index1", "")), len(sample.get("Index2", "")))))
            # List of sub-demux with a decreasing order of index lengths
            sub_demux_list = [x[0] for x in sorted(sub_demux_index_lens, key=lambda x: sum(x[1]), reverse=True)]
            sub_demux_with_max_index_lens = sub_demux_list[0]
            # Start with the unassigned list with the longest index
            max_unassigned_csv = os.path.join(self.run_dir, f"Demultiplexing_{sub_demux_with_max_index_lens}", "UnassignedSequences.csv")
            with open(max_unassigned_csv, 'r') as max_unassigned_file:
                reader = csv.DictReader(max_unassigned_file)
                max_unassigned_indexes = [row for row in reader]
            # Filter by lane
            max_unassigned_indexes = [idx for idx in max_unassigned_indexes if idx["Lane"] == lane]
            # Complicated case with multiple demuxes. Take the full list if there is only one sub-demux otherwise
            if len(sub_demux_list) > 1:
                # Order: from longer to shorter indexes
                sub_demux_with_shorter_index_lens = sub_demux_list[1:]
                for sub_demux in sub_demux_with_shorter_index_lens:
                    unassigned_csv = os.path.join(self.run_dir, f"Demultiplexing_{sub_demux}", "UnassignedSequences.csv")
                    with open(unassigned_csv, 'r') as unassigned_file:
                        reader = csv.DictReader(unassigned_file)
                        unassigned_indexes = [row for row in reader]
                    # Filter by lane
                    unassigned_indexes = [unassigned_index for unassigned_index in unassigned_indexes if unassigned_index["Lane"] == lane]
                    # Remove overlapped indexes from the list of max_unassigned_indexes
                    idx1_overlapped_len = min([demux_lens_pair[1] for demux_lens_pair in sub_demux_index_lens if demux_lens_pair[0] == sub_demux][0][0],
                                              [demux_lens_pair[1] for demux_lens_pair in sub_demux_index_lens if demux_lens_pair[0] == sub_demux_with_max_index_lens][0][0])
                    idx2_overlapped_len = min([demux_lens_pair[1] for demux_lens_pair in sub_demux_index_lens if demux_lens_pair[0] == sub_demux][0][1],
                                              [demux_lens_pair[1] for demux_lens_pair in sub_demux_index_lens if demux_lens_pair[0] == sub_demux_with_max_index_lens][0][1])
                    for unassigned_index in unassigned_indexes:
                        idx1_overlapped_seq = unassigned_index['I1'][:idx1_overlapped_len]
                        idx2_overlapped_seq = unassigned_index['I2'][:idx2_overlapped_len]
                        # Remove the overlapped record from the max_unassigned_indexes list
                        max_unassigned_indexes = [max_unassigned_index for max_unassigned_index in max_unassigned_indexes if not (max_unassigned_index['I1'][:idx1_overlapped_len] == idx1_overlapped_seq and max_unassigned_index['I2'][:idx2_overlapped_len] == idx2_overlapped_seq)]
            # Append to the aggregated_unassigned_indexes list
            aggregated_unassigned_indexes += max_unassigned_indexes
        # Sort aggregated_unassigned_indexes list first by lane and then by Count in the decreasing order
        aggregated_unassigned_indexes = sorted(aggregated_unassigned_indexes, key=lambda x: (x['Lane'], -int(x['Count'])))
        # Write to a new UnassignedSequences.csv file under
        aggregated_unassigned_csv = os.path.join(self.run_dir, self.demux_dir, "UnassignedSequences.csv")
        write_to_csv(aggregated_unassigned_indexes, aggregated_unassigned_csv)


    # Aggregate demux results
    def aggregate_demux_results(self, demux_results_dirs):
        # In case of single demux
        if len(demux_results_dirs) == 1:
            # TODO: Check NoIndex case. Can Base2Fastq generate FastQs for both reads and indexes for NoIndex sample?
            # Otherwise just softlink contents of Demultplexing_0 into Demultiplexing
            symlink_demux_dir(demux_results_dirs[0], os.path.join(self.run_dir, self.demux_dir))
        else:
            # Ensure the destination directory exists
            if not os.path.exists(os.path.join(self.run_dir, self.demux_dir):
                os.makedirs(os.path.join(self.run_dir, self.demux_dir)
            # Clear all content under dest_dir
            clear_dir(os.path.join(self.run_dir, self.demux_dir)
            demux_runmanifest = collect_demux_runmanifest(demux_results_dirs)
            # Aggregate the output FastQ files of samples from multiple demux
            aggregate_sample_fastq(demux_runmanifest)
            # Symplink the output FastQ files of undet only if a lane does not have multiple demux
            aggregate_undet_fastq(demux_runmanifest)
            # Aggregate stats in IndexAssignment.csv
            TBD
            # Aggregate stats in UnassignedSequences.csv
            aggregate_stats_unassigned(demux_runmanifest)
            # Aggregate stats in Project_RunStats.json
            TBD



    def upload_demux_results_to_statusdb(self):
        # TODO: dump contents of IndexAssignment.csv and UnassignedSequences.csv into statusdb document
        doc_obj = self.db.get_db_entry(self.NGI_run_id)
        index_assignement_file = os.path.join(self.run_dir, "Demultiplexing", "IndexAssignment.csv")
        with open(index_assignement_file, 'r') as index_file:
            reader = csv.DictReader(index_file)
            index_assignments = [row for row in reader]
        unassigned_sequences_file = os.path.join(self.run_dir, "Demultiplexing", "UnassignedSequences.csv")
        with open(unassigned_sequences_file, 'r') as unassigned_file:
            reader = csv.DictReader(unassigned_file)
            unassigned_sequences = [row for row in reader]
        project_dirs = [f.path for f in os.scandir(os.path.join(self.run_dir, "Demultiplexing")) if f.is_dir() and not "PhiX" in f]
        for project_dir in project_dirs:
            run_stats_file = glob.glob(os.path.join(project_dir, "*_RunStats.json"))
            with open(run_stats_file) as stats_json:
                project_sample_stats_raw = json.load(stats_json)
            collected_sample_stats = {}
            for sample_stats in project_sample_stats_raw["SampleStats"]:
                sample_name = sample_stats["SampleName"]
                percent_q30 = sample_stats["PercentQ30"]
                quality_score_mean = sample_stats["QualityScoreMean"]
                percent_mismatch = sample_stats["PercentMismatch"]
                collected_sample_stats[sample_name] = {
                    "PercentQ30": percent_q30,
                    "QualityScoreMean": quality_score_mean,
                    "PercentMismatch": percent_mismatch
                    }
            for assignment in index_assignments:
                sample = assignment.get("SampleName")
                sample_stats_to_add = collected_sample_stats.get(sample)
                assignment["PercentQ30"] = sample_stats_to_add.get("PercentQ30")
                assignment["QualityScoreMean"] = sample_stats_to_add.get("QualityScoreMean")
                assignment["PercentMismatch"] = sample_stats_to_add.get("PercentMismatch")
        demultiplex_stats = {
            "Demultiplex_Stats": {
                "Index_Assignment": index_assignments,
                "Unassigned_Sequences": unassigned_sequences
                }
            }
        doc_obj["Aviti": demultiplex_stats]
        self.db.upload_to_statusdb(doc_obj)

    def sync_metadata(self):
        # TODO: copy metadata from demuxed run to ngi-nas-ns
        pass

    def make_transfer_indicator(self):
        transfer_indicator = os.path.join(self.run_dir, '.rsync_ongoing')
        Path(transfer_indicator).touch()

    def transfer(self):
        transfer_details = self.CONFIG.get("Element").get(self.sequencer_type).get("transfer_details") #TODO: Add section to taca.yaml
        command = ("rsync"
                   + " -rLav"
                   + f" --chown={transfer_details.get("owner")}"
                   + f" --chmod={transfer_details.get("permissions")}"
                   + " --exclude BaseCalls" # TODO: check that we actually want to exclude these
                   + " --exclude Alignment"
                   + f" {self.run_dir}"
                   + f" {transfer_details.get("user")@transfer_details.get("host")}:/"
                   + "; echo $? > .rsync_exit_status"
            )  # TODO: any other options?
        try:
            p_handle = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            logger.info(
                "Transfer to analysis cluster "
                f"started for run {self} on {datetime.now()}"
            )
        except subprocess.CalledProcessError:
            logger.warning("An error occurred while starting transfer to analysis cluster "
                            f"for {self} on {datetime.now()}."
            )
        return

    def remove_transfer_indicator(self):
        # TODO: remove hidden file in run directory
        pass

    def update_transfer_log(self):
        # TODO: update the transfer log
        pass

    def archive(self):
        # TODO: move run dir to nosync
        pass
