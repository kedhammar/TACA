from taca.nanopore import instrument_transfer
from unittest.mock import patch, mock_open, call, Mock
import tempfile
import pytest
import os

DUMMY_RUN_NAME = "20240112_2342_MN19414_TEST12345_randomhash"

# To check coverage, use
# pytest --cov=taca.nanopore.instrument_transfer test_instrument_transfer.py

def test_sequencing_finished():
    with patch("os.listdir") as mock_listdir:
        mock_listdir.return_value = ["file1", "file2", "final_summary"]
        assert instrument_transfer.sequencing_finished("path") is True

        mock_listdir.return_value = ["file1", "file2"]
        assert instrument_transfer.sequencing_finished("path") is False


def test_dump_path():
    with patch("builtins.open", new_callable=mock_open) as mock_file:
        instrument_transfer.dump_path("path/to/run")
        mock_file.assert_called_once_with("path/to/run/run_path.txt", "w")


def test_write_finished_indicator():
    with patch("builtins.open", new_callable=mock_open) as mock_file:
        run_path = "path/to/run"
        instrument_transfer.write_finished_indicator(run_path)
        mock_file.assert_called_once_with(run_path + "/.sync_finished", "w")


def test_sync_to_storage():
    with patch("subprocess.Popen") as mock_Popen:
        instrument_transfer.sync_to_storage("run_dir", "destination", "log")
        mock_Popen.assert_called_once_with(
            [
                "run-one",
                "rsync",
                "-rvu",
                "--log-file=" + "log",
                "run_dir",
                "destination",
            ]
        )


@patch("taca.nanopore.instrument_transfer.archive_finished_run")
@patch("taca.nanopore.instrument_transfer.write_finished_indicator")
@patch("subprocess.run")
def test_final_sync_to_storage(
    mock_run,
    mock_write_finished_indicator,
    mock_archive_finished_run,
):
    mock_write_finished_indicator.return_value = ".sync_finished"

    # For finished run
    mock_run.return_value.returncode = 0

    instrument_transfer.final_sync_to_storage(
        run_dir="run_dir",
        destination="destination",
        archive_dir="archive_dir",
        log="log_path",
    )

    assert mock_run.call_args_list[0] == call(
        [
            "run-one",
            "rsync",
            "-rvu",
            "--log-file=" + "log_path",
            "run_dir",
            "destination",
        ]
    )

    assert mock_run.call_args_list[1] == call(
        ["rsync", ".sync_finished", "destination/run_dir"]
    )

    mock_archive_finished_run.assert_called_once_with("run_dir", "archive_dir")

    # For not finished run
    mock_run.return_value.returncode = 1

    instrument_transfer.final_sync_to_storage(
        run_dir="run_dir",
        destination="destination",
        archive_dir="archive_dir",
        log="log_path",
    )

    assert mock_run.call_count == 3


def test_archive_finished_run():
    # Set up combinatorial testing for cases
    archive_dirs = [
        "/data/nosync",
        "/data/nosync/experiment/sample",
        "/data/nosync/experiment/sample/run",
    ]
    neighbor_dirs = [
        None,
        "/data/experiment/sample/run2",
        "/data/experiment/sample2",
    ]

    for archive_dir in archive_dirs:
        for neighbor_dir in neighbor_dirs:
            # Set up tmp dir
            tmp = tempfile.TemporaryDirectory()
            tmp_path = tmp.name

            # Create run dir
            experiment_path = tmp_path + "/data/experiment"
            sample_path = experiment_path + "/sample"
            run_path = sample_path + f"/{DUMMY_RUN_NAME}"
            os.makedirs(run_path)

            # Create neightbor dir, if any
            if neighbor_dir:
                neighbor_path = tmp_path + neighbor_dir
                os.mkdir(neighbor_path)

            # Create archive dir
            archive_path = tmp_path + archive_dir
            os.makedirs(archive_path)

            # Execute code
            instrument_transfer.archive_finished_run(run_path, archive_path)

            # Assert run is moved to archive dir
            assert os.path.exists(archive_path + f"/experiment/sample/{DUMMY_RUN_NAME}")

            # Assert run is removed from original location
            assert not os.path.exists(run_path)

            # Assert experiment and sample dirs are removed if empty
            if neighbor_dir:
                assert os.path.exists(neighbor_path)
                if neighbor_dir == "/data/experiment/sample/run2":
                    assert os.path.exists(sample_path)
                else:
                    assert not os.path.exists(sample_path)
            else:
                assert not os.path.exists(experiment_path)

            tmp.cleanup()


def create_logs() -> list:
    tmp = tempfile.TemporaryDirectory()

    # For each instrument position dir and corresponding flowcell name...
    for position_dir, flowcell_name in {
        "1A": "PAM12345",
        "MN19414": "FLG12345",
    }.items():
        # --> Create position dir
        position_path = tmp.name + f"/log/{position_dir}"
        os.makedirs(position_path)

        # --> Populate each position dir with two log files
        for file_n in range(1, 3):
            # --> Populate each log file with two entries
            lines = []
            for entry_n in range(1, 3):
                lines += [
                    "2023-10-31 15:12:54.1354    INFO: mux_scan_result (user_messages)",
                    f"    flow_cell_id: {flowcell_name}",
                    f"    num_pores: {file_n}",
                    f"    total_pores: {entry_n}",
                ]

            with open(position_path + f"/control_server_log-{file_n}.txt", "w") as file:
                file.write("\n".join(lines))

    # Run code
    logs = instrument_transfer.parse_position_logs(tmp.name + "/log")

    return logs


def test_parse_position_logs():
    logs = create_logs()

    # Create template to compare output against
    template = []
    for position_dir, flowcell_name in {
        "MN19414": "FLG12345",
        "1A": "PAM12345",
    }.items():
        for file_n in range(1, 3):
            template += [
                {
                    "position": position_dir,
                    "timestamp": "2023-10-31 15:12:54.1354",
                    "category": "INFO: mux_scan_result (user_messages)",
                    "body": {
                        "flow_cell_id": flowcell_name,
                        "num_pores": str(file_n),
                        "total_pores": str(entry_n),
                    },
                }
                for entry_n in range(1, 3)
            ]

    assert logs == template


def test_get_pore_counts():
    logs = create_logs()
    pore_counts = instrument_transfer.get_pore_counts(logs)

    # Create template to compare output against
    template = []
    for pos, fc in {
        "MN19414": "FLG12345",
        "1A": "PAM12345",
    }.items():
        for i in range(1, 3):
            for j in range(1, 3):
                template.append(
                    {
                        "flow_cell_id": fc,
                        "timestamp": "2023-10-31 15:12:54.1354",
                        "position": pos,
                        "type": "mux",
                        "num_pores": str(i),
                        "total_pores": str(j),
                    }
                )

    assert pore_counts == template


def test_dump_pore_count_history():
    logs = create_logs()
    pore_counts = instrument_transfer.get_pore_counts(logs)

    # Nothing to add, no file
    tmp = tempfile.TemporaryDirectory()
    run_path = tmp.name + f"/experiment/sample/{DUMMY_RUN_NAME}"
    os.makedirs(run_path)
    new_file = instrument_transfer.dump_pore_count_history(run_path, pore_counts)
    assert open(new_file, "r").read() == ""
    tmp.cleanup()

    # Nothing to add, file is present
    tmp = tempfile.TemporaryDirectory()
    run_path = tmp.name + f"/experiment/sample/{DUMMY_RUN_NAME}"
    os.makedirs(run_path)
    open(run_path + "/pore_count_history.csv", "w").write("test")
    new_file = instrument_transfer.dump_pore_count_history(run_path, pore_counts)
    assert open(new_file, "r").read() == "test"
    tmp.cleanup()

    # Something to add
    tmp = tempfile.TemporaryDirectory()
    run_path = tmp.name + f"/experiment/sample/{DUMMY_RUN_NAME.replace('TEST','FLG')}"
    os.makedirs(run_path)
    new_file = instrument_transfer.dump_pore_count_history(run_path, pore_counts)

    template = (
        "\n".join(
            [
                "flow_cell_id,timestamp,position,type,num_pores,total_pores",
                "FLG12345,2023-10-31 15:12:54.1354,MN19414,mux,1,1",
                "FLG12345,2023-10-31 15:12:54.1354,MN19414,mux,1,2",
                "FLG12345,2023-10-31 15:12:54.1354,MN19414,mux,2,1",
                "FLG12345,2023-10-31 15:12:54.1354,MN19414,mux,2,2",
            ]
        )
        + "\n"
    )

    assert open(new_file, "r").read() == template
    tmp.cleanup()
