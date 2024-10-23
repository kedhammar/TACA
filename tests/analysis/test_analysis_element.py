import logging
from io import StringIO
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd
import pytest
from dirhash import dirhash

from tests.element.test_Element_Runs import create_element_run_dir, get_config


def parametrize_testruns() -> tuple[list[dict], list[str]]:
    """Helper function to build test parametrization from a friendly string table."""

    testrun_descs: list[str] = ["ready to demux", "demux_ongoing"]

    kwarg_table = """
    lims_manifest  metadata_files  run_finished  outcome_completed  demux_dir  demux_done  rsync_ongoing  rsync_exit_status  nosync
    True           True            True          True               False      False       False          None               False
    True           True            True          True               True       False       False          None               False
    """

    # Turn string table to datastream
    data = StringIO(kwarg_table)

    # Read data, trimming whitespace
    df = pd.read_csv(data, sep=r"\s+")

    # Compile into list of parameters to use
    testrun_kwargs: list[dict] = df.to_dict(orient="records")

    assert len(testrun_descs) == len(testrun_kwargs)

    return testrun_kwargs, testrun_descs


testrun_kwargs, testrun_descs = parametrize_testruns()


@pytest.mark.parametrize("run_kwargs", testrun_kwargs, ids=testrun_descs)
def test_run_preprocessing(create_dirs, run_kwargs):
    tmp: TemporaryDirectory = create_dirs

    # Mock config
    config = get_config(tmp)
    mock_config = patch("taca.utils.config.CONFIG", new=config)
    mock_config.start()

    # Mock DB
    mock_db = patch("taca.element.Element_Runs.ElementRunsConnection")
    mock_db.start()

    # Mock subprocess
    mock_subprocess = patch("subprocess.Popen")
    mock_subprocess.start()

    # Create run dir and associated LIMS manifest
    run_dir = create_element_run_dir(tmp=tmp, **run_kwargs)

    # Import module to test
    from taca.analysis import analysis_element as to_test

    # Test
    to_test.run_preprocessing(run_dir)

    # Stop mocks
    patch.stopall()


@pytest.fixture
def aviti_fixture(create_dirs, caplog):
    # Create tempdir
    tmp: TemporaryDirectory = create_dirs

    # Capture log
    caplog.at_level(logging.INFO)

    # Mocks
    mocks = {
        "mock_config": patch("taca.utils.config.CONFIG", new=get_config(tmp)).start(),
        "mock_db": patch(
            "taca.element.Element_Runs.ElementRunsConnection", autospec=True
        ).start(),
        "mock_mail": patch("taca.analysis.analysis_element.send_mail").start(),
        "mock_popen": patch("subprocess.Popen").start(),
    }

    # Import module to test
    from taca.analysis import analysis_element as to_test

    # Yield fixtures
    yield to_test, tmp, caplog, mocks

    # Stop mocks
    patch.stopall()


def test_process_on_empty_dir(aviti_fixture):
    to_test, tmp, caplog, mocks = aviti_fixture
    """Should raise FileNotFoundError when no files are present in the run dir and send mail."""

    # Create dir
    run_dir = create_element_run_dir(
        tmp=tmp,
        lims_manifest=False,
        metadata_files=False,
        run_finished=False,
        outcome_completed=False,
        demux_dir=False,
        demux_done=False,
        rsync_ongoing=False,
        rsync_exit_status=None,
        nosync=False,
    )

    with pytest.raises(FileNotFoundError):
        to_test.run_preprocessing(run_dir)

    # Assertions
    mocks["mock_mail"].assert_called_once()
    assert "Run parameters file not found" in caplog.text


def test_process_on_dir_w_metadata(aviti_fixture):
    """Should update statusdb."""
    to_test, tmp, caplog, mocks = aviti_fixture

    # Sub-mock configuration
    mocks["mock_db"].return_value.check_db_run_status.return_value = "ongoing"
    mocks["mock_db"].return_value.upload_to_statusdb.return_value = None

    # Add metadata files
    run_dir = create_element_run_dir(
        tmp=tmp,
        lims_manifest=False,
        metadata_files=True,
        run_finished=False,
        outcome_completed=False,
        demux_dir=False,
        demux_done=False,
        rsync_ongoing=False,
        rsync_exit_status=None,
        nosync=False,
    )

    to_test.run_preprocessing(run_dir)

    assert mocks["mock_db"].return_value.upload_to_statusdb.called


@pytest.mark.skip("Currently a failed run is treated as an ongoing run.")
def test_process_on_failed_run(aviti_fixture):
    to_test, tmp, caplog, mocks = aviti_fixture

    # Sub-mock configuration
    mocks["mock_db"].return_value.check_db_run_status.return_value = "ongoing"
    mocks["mock_db"].return_value.upload_to_statusdb.return_value = None

    # Add metadata files
    run_dir = create_element_run_dir(
        tmp=tmp,
        lims_manifest=False,
        metadata_files=True,
        run_finished=True,
        outcome_completed=False,
        demux_dir=False,
        demux_done=False,
        rsync_ongoing=False,
        rsync_exit_status=None,
        nosync=False,
    )

    to_test.run_preprocessing(run_dir)


def test_process_on_finished_run_wo_lims_manifest(aviti_fixture):
    """Should fail to find LIMS run manifest and send mail."""
    to_test, tmp, caplog, mocks = aviti_fixture

    # Sub-mock configuration
    mocks["mock_db"].return_value.check_db_run_status.return_value = "ongoing"
    mocks["mock_db"].return_value.upload_to_statusdb.return_value = None

    # Add metadata files
    run_dir = create_element_run_dir(
        tmp=tmp,
        lims_manifest=False,
        metadata_files=True,
        run_finished=True,
        outcome_completed=True,
        demux_dir=False,
        demux_done=False,
        rsync_ongoing=False,
        rsync_exit_status=None,
        nosync=False,
    )

    to_test.run_preprocessing(run_dir)

    assert "No manifest found for run" in caplog.text
    mocks["mock_mail"].assert_called_once()


def test_process_on_finished_run(aviti_fixture):
    """Should start demux."""
    to_test, tmp, caplog, mocks = aviti_fixture

    # Sub-mock configuration
    mocks["mock_db"].return_value.check_db_run_status.return_value = "ongoing"
    mocks["mock_db"].return_value.upload_to_statusdb.return_value = None

    # Add metadata files
    run_dir = create_element_run_dir(
        tmp=tmp,
        lims_manifest=True,
        metadata_files=True,
        run_finished=True,
        outcome_completed=True,
        demux_dir=False,
        demux_done=False,
        rsync_ongoing=False,
        rsync_exit_status=None,
        nosync=False,
    )

    to_test.run_preprocessing(run_dir)

    expected_call = " ".join(
        [
            "mock_bases2fastq_path",
            f"{tmp.name}/ngi_data/sequencing/AV242106/20240926_AV242106_A2349523513",
            f"{tmp.name}/ngi_data/sequencing/AV242106/20240926_AV242106_A2349523513/Demultiplexing_0",
            "-p 8",
            "--num-unassigned 500",
            f"-r {tmp.name}/ngi_data/sequencing/AV242106/20240926_AV242106_A2349523513/20240926_AV242106_A2349523513_demux_0.csv",
            "--legacy-fastq",
            "--force-index-orientation",
        ]
    )
    assert mocks["mock_popen"].call_args.args[0] == expected_call
    assert "Bases2Fastq conversion and demultiplexing started for run " in caplog.text
    assert mocks["mock_db"].return_value.upload_to_statusdb.called
