from taca.nanopore import instrument_transfer
from unittest.mock import patch, mock_open
import tempfile


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


# TODO WIP
def todo_test_archive_finished_run():

    instrument_transfer.archive_finished_run("/data/experiment/sample/run", "/data/nosync")
    # Creating directory
    # 1. Archive is empty --> Create experiment dir and sample dir
    # 2. Archive has matching experiment dir --> Create sample dir
    # 3. Archive has matching experiment dir and sample dir --> Do nothing

    # Cleaning up
    # 1. Sample folder is not empty --> Keep it
    # 2. Sample folder is empty --> Remove it
    # 3. Experiment folder is not empty --> Keep it
    # 4. Experiment folder is empty --> Remove it