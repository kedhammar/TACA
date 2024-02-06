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

    # Set up config mocks
    patch("taca.utils.config.CONFIG", new=test_config_yaml).start()
    patch("taca.nanopore.ONT_run_classes.CONFIG", new=test_config_yaml).start()

    # Set up mocks
    mock_db = patch(
        "taca.nanopore.ONT_run_classes.NanoporeRunsConnection",
    ).start()
    mock_db.return_value.check_run_exists.return_value = False
    mock_db.return_value.check_run_status.return_value = "ongoing"

    mock_error = patch("taca.analysis.analysis_nanopore.send_error_mail").start()

    # Reload module to add mocks
    importlib.reload(analysis_nanopore)

    # Create run dir
    create_run_dir(tmp, script_files=True, run_finished=True, sync_finished=True)

    # Start testing
    analysis_nanopore.ont_transfer(run_abspath=None, qc=False)
    print(mock_error.call_args_list)

    for record in caplog.records:
        print(record.levelname, record.message)
