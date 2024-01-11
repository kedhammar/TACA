from taca.nanopore import instrument_transfer
from unittest.mock import patch, mock_open, call, Mock
import tempfile
import pytest


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

@pytest.mark.skip
def test_archive_finished_run():
    # TODO
    pass
    # Creating directory
    # 1. Archive is empty --> Create experiment dir and sample dir
    # 2. Archive has matching experiment dir --> Create sample dir
    # 3. Archive has matching experiment dir and sample dir --> Do nothing

    # Cleaning up
    # 1. Sample folder is not empty --> Keep it
    # 2. Sample folder is empty --> Remove it
    # 3. Experiment folder is not empty --> Keep it
    # 4. Experiment folder is empty --> Remove it

@pytest.mark.skip
def test_parse_position_logs():
    # TODO
    pass

@pytest.mark.skip
def test_get_pore_counts():
    # TODO
    pass

@pytest.mark.skip
def test_dump_pore_count_history():
    # TODO
    pass