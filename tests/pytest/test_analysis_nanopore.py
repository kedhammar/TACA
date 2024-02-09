import importlib
import subprocess
from unittest.mock import patch

from test_ONT_run_classes import (
    create_run_dir,
    make_test_config,
)

from taca.analysis import analysis_nanopore


def test_ont_transfer(create_dirs, caplog):
    """Test the "taca analaysis ont-transfer" subcommand automation from
    start to finish for a variety of runs.
    """

    # Create dir tree from fixture
    tmp = create_dirs

    # Create test config
    test_config_yaml = make_test_config(tmp)

    ## MOCKS

    # Mock config
    patch("taca.utils.config.CONFIG", new=test_config_yaml).start()
    patch("taca.nanopore.ONT_run_classes.CONFIG", new=test_config_yaml).start()

    # Mock database connection
    mock_db = patch(
        "taca.nanopore.ONT_run_classes.NanoporeRunsConnection",
    ).start()
    mock_db.return_value.check_run_exists.return_value = False
    mock_db.return_value.check_run_status.return_value = "ongoing"
    mock_db.return_value.finish_ongoing_run

    # Mock parsing MinKNOW auxillary files
    patch("taca.nanopore.ONT_run_classes.ONT_run.parse_minknow_json").start()
    patch("taca.nanopore.ONT_run_classes.ONT_run.parse_pore_activity").start()

    # Mock subprocess.Popen ONLY for Anglerfish
    original_popen = subprocess.Popen

    def side_effect(*args, **kwargs):
        if "anglerfish" in args[0]:
            return mock_Popen
        else:
            return original_popen(*args, **kwargs)

    mock_Popen = patch(
        "taca.nanopore.ONT_run_classes.subprocess.Popen", side_effect=side_effect
    ).start()
    mock_Popen.pid = 1337  # Nice

    # Reload module to implement mocks
    importlib.reload(analysis_nanopore)

    ## CREATE RUN DIRS

    # User run
    create_run_dir(
        tmp,
        run_id="TestUserRun",
        script_files=True,
        run_finished=True,
        sync_finished=True,
    )
    # QC run
    create_run_dir(
        tmp,
        qc=True,
        run_id="TestQCRun",
        script_files=True,
        run_finished=True,
        sync_finished=True,
        anglerfish_samplesheets=True,
        fastq_dirs=True,
        barcode_dirs=True,
    )

    # Start testing
    analysis_nanopore.ont_transfer(run_abspath=None, qc=False)
