""" Transfers new PromethION runs to ngi-nas using rsync.
"""
__version__ = "1.0.5"

import os
import re
import shutil
import pathlib
import argparse
import subprocess
from glob import glob

def main(args):
    """Find promethion runs and transfer them to storage. 
    Archives the run when the transfer is complete."""
    data_dir = args.source_dir
    run_pattern = re.compile("\d{8}_\d{4}_[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+")
    destination_dir = args.dest_dir
    archive_dir = args.archive_dir
    log_file = os.path.join(data_dir, 'rsync_log.txt')

    position_logs = parse_position_logs(args.minknow_logs_dir)
    qc_logs = get_qc_logs(position_logs)

    runs = [
        path
        for path in glob("*/*/*", recursive=True)
        if re.match(run_pattern, os.path.basename(path))
    ]
    
    # Split finished and unfinished runs
    not_finished = []
    finished = []
    
    for run in runs:
        if sequencing_finished(run):
            finished.append(run)
        else:
            not_finished.append(run)
        dump_path(run)
        dump_qc_logs(run, qc_logs)

    # Start transfer of unfinished runs first (detatched)
    for run in not_finished:
        sync_to_storage(run, destination_dir, log_file)
    for run in finished:
        final_sync_to_storage(run, destination_dir, archive_dir, log_file) 

def sequencing_finished(run_dir):
    sequencing_finished_indicator = 'final_summary'
    run_dir_content = os.listdir(run_dir) 
    for item in run_dir_content:
        if sequencing_finished_indicator in item:
            return True
    return False

def dump_path(run_path):
    """Dump path to run to a file that can be
    used when uploading stats to statusdb from preproc."""
    new_file = os.path.join(run_path, 'run_path.txt')
    proj, sample, run = run_path.split('/')[-3:]
    path_to_write = os.path.join(proj, sample, run)
    with open(new_file, 'w') as f:
        f.write(path_to_write)
    
def write_finished_indicator(run_path):
    """Write a hidden file to indicate 
    when the finial rsync is finished."""
    new_file = os.path.join(run_path, '.sync_finished')
    pathlib.Path(new_file).touch()
    return new_file

def sync_to_storage(run_dir, destination, log_file):
    """Sync the run to storage using rsync. 
    Skip if rsync is already running on the run."""
    command = ['run-one', 'rsync', '-rv', '--log-file=' + log_file, run_dir, destination]
    process_handle = subprocess.Popen(command)
    print('Initiated rsync with the following parameters: {}'.format(command))
    
def final_sync_to_storage(run_dir, destination, archive_dir, log_file):
    """Do a final sync of the run to storage, then archive it. 
    Skip if rsync is already running on the run."""
    print('Performing a final sync of {} to storage'.format(run_dir))
    command = ['run-one', 'rsync', '-rv', '--log-file=' + log_file, run_dir, destination]
    process_handle = subprocess.run(command)
    if process_handle.returncode == 0:
        finished_indicator = write_finished_indicator(run_dir)
        dest = os.path.join(destination, os.path.basename(run_dir))
        sync_finished_indicator = ['rsync', finished_indicator, dest]
        process_handle = subprocess.run(sync_finished_indicator)
        archive_finished_run(run_dir, archive_dir)
    else:
        print('Previous rsync might be running still. Skipping {} for now.'.format(run_dir))
        return

def archive_finished_run(run_dir, archive_dir):
    """Move finished run to archive (nosync)."""
    dir_to_move = str(pathlib.Path(run_dir).parent)
    print('Archiving {}'.format(dir_to_move))
    top_dir = str(pathlib.Path(run_dir).parent.parent)
    project_id = os.path.basename(top_dir)
    project_archive = os.path.join(archive_dir, project_id)
    if os.path.exists(project_archive):
        shutil.move(dir_to_move, project_archive)
    else:
        os.mkdir(project_archive)
        shutil.move(dir_to_move, project_archive)
    if not os.listdir(top_dir):
        print("Project folder {} is empty. Removing it.".format(top_dir))
        os.rmdir(top_dir)
    else:
        print("Some data is still left in {}. Keeping it.".format(top_dir))  # Might be another run for the same project


def parse_position_logs(minknow_log_dir: str) -> list:
    """Look through all flow cell position logs and boil down into a strucutred list of dicts

    Example output:
    [{
        "timestamp": "2023-07-10 15:44:31.481512",
        "category": "INFO: platform_qc.report (user_messages)",
        "body": {
            "flow_cell_id": "PAO33763"
            "num_pores": "8378"
        }
    } ... ]

    """

    position_log_file_name = "control_server_log-0.txt"
    log_timestamp_format = "%Y-%m-%d %H:%M:%S.%f"

    positions = []
    for col in "123":
        for row in "ABCDEFGH":
            positions.append(col + row)

    entries = []
    for position in positions:
        with open(
            f"{minknow_log_dir}/{position}/{position_log_file_name}", "r"
        ) as stream:
            lines = stream.readlines()
            for i in range(0, len(lines)):
                line = lines[i]
                if line[0:4] != "    ":
                    # Line is log header
                    split_header = line.split(" ")
                    timestamp = " ".join(split_header[0:2])
                    category = " ".join(split_header[2:])

                    entry = {
                        "position": position,
                        "timestamp": timestamp.strip(),
                        "category": category.strip(),
                    }
                    entries.append(entry)
                else:
                    # Line is log body
                    if "body" not in entry:
                        entry["body"] = {}
                    key = line.split(": ")[0].strip()
                    val = ": ".join(line.split(": ")[1:]).strip()
                    entry["body"][key] = val

    return entries


def get_qc_logs(position_logs: list) -> list:
    """Take the flowcell log list output by parse_position_logs() and subset to contain only QC info."""

    qc_logs = []
    for entry in position_logs:
        if "INFO: platform_qc.report (user_messages)" in entry["category"]:

            new_entry = {
                "flow_cell_id": entry["body"]["flow_cell_id"],
                "timestamp": entry["timestamp"],
                "position": entry["position"],
                "num_pores": entry["body"]["num_pores"],
            }

            qc_logs.append(new_entry)

    return qc_logs


def dump_qc_logs(run_path, qc_logs):

    flowcell_id = os.path.basename(run_path).split["_"][-2]
    new_file_path = os.path.join(run_path, "flowcell_qc_logs.csv")

    flowcell_qc_logs = [
        log_entry for log_entry in qc_logs if log_entry["flow_cell_id"] == flowcell_id
    ]

    flowcell_qc_logs_sorted = sorted(
        flowcell_qc_logs,
        key=lambda x: (
            x["flow_cell_id"],
            x["timestamp"],
            x["position"],
            x["num_pores"],
        ),
    )

    header = flowcell_qc_logs_sorted[0].keys()
    rows = [e.values() for e in flowcell_qc_logs_sorted]

    with open(new_file_path, "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(row) + "\n")


if __name__ == "__main__":
    # This is clunky but should be fine since it will only ever run as a cronjob
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        dest="source_dir",
        help="Full path to directory containing runs to be synced.",
    )
    parser.add_argument(
        "--dest",
        dest="dest_dir",
        help="Full path to destination directory to sync runs to.",
    )
    parser.add_argument(
        "--archive",
        dest="archive_dir",
        help="Full path to directory containing runs to be synced.",
    )
    parser.add_argument(
        "--minknow_logs",
        dest="minknow_logs_dir",
        help="Full path to the directory containing the MinKNOW position logs.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()
    
    main(args)