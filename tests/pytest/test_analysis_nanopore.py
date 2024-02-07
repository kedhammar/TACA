import importlib
from unittest.mock import patch

from test_ONT_run_classes import (
    create_run_dir,
    make_test_config,
)

from taca.analysis import analysis_nanopore

# To check coverage, use
# pytest -s --log-cli-level=INFO --cov=taca.analysis.analysis_nanopore --cov-report term-missing -vv tests/pytest/test_analysis_nanopore.py


def test_ont_transfer(create_dirs, caplog):
    """Test the "taca analaysis ont-transfer" subcommand automation from
    start to finish for a variety of runs.
    """

    # Create dir tree from fixture
    tmp = create_dirs

    # Create test config
    test_config_yaml = make_test_config(tmp)

    # === MOCKS =============================================================

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
    patch("taca.nanopore.ONT_run_classes.ONT_user_run.parse_minknow_json").start()
    patch("taca.nanopore.ONT_run_classes.ONT_user_run.parse_pore_activity").start()

    # Reload module to implement mocks
    importlib.reload(analysis_nanopore)

    # === CREATE RUN DIRS ===================================================

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
    )

    # Start testing
    analysis_nanopore.ont_transfer(run_abspath=None, qc=False)
