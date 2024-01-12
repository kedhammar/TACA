from taca.nanopore import instrument_transfer
from unittest.mock import patch, mock_open, call, Mock
import tempfile
import pytest
import os

DUMMY_RUN_NAME = "20240112_2342_MN19414_TEST12345_randomhash"


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

            # Assert experiment and sample dirs are removed if appropriate
            if neighbor_dir:
                if neighbor_dir == "/data/experiment/sample/run2":
                    assert os.path.exists(sample_path)
                else:
                    assert not os.path.exists(sample_path)
            else:
                assert not os.path.exists(experiment_path)

            tmp.cleanup()


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
