import importlib
import logging
import subprocess
from io import StringIO
from unittest.mock import patch

import pandas as pd
import pytest

from taca.analysis import analysis_nanopore
from tests.nanopore.test_ONT_run_classes import (
    create_ONT_run_dir,
    make_ONT_test_config,
)


def build_run_properties() -> dict:
    """In order to parametrize the test in a comprehensive way, the parametrization is
    tabulated as a string here.
    """

    parameter_string_table = """
    instrument qc    run_finished sync_finished raw_dirs fastq_dirs barcode_dirs anglerfish_samplesheets anglerfish_ongoing anglerfish_exit
    promethion False False        False         False    False      False        False                   False              None
    promethion False True         False         False    False      False        False                   False              None
    promethion False True         True          False    False      False        False                   False              None
    promethion False True         True          True     False      False        False                   False              None
    promethion False True         True          True     True       False        False                   False              None
    promethion False True         True          True     True       True         False                   False              None
    minion     False False        False         False    False      False        False                   False              None
    minion     False True         False         False    False      False        False                   False              None
    minion     False True         True          False    False      False        False                   False              None
    minion     False True         True          True     False      False        False                   False              None
    minion     False True         True          True     True       False        False                   False              None
    minion     False True         True          True     True       True         False                   False              None
    minion     True  False        False         False    False      False        False                   False              None
    minion     True  True         False         False    False      False        False                   False              None
    minion     True  True         True          False    False      False        False                   False              None
    minion     True  True         True          True     False      False        False                   False              None
    minion     True  True         True          True     True       False        False                   False              None
    minion     True  True         True          True     True       True         False                   False              None
    minion     True  True         True          True     True       True         True                    False              None
    minion     True  True         True          True     True       True         True                    True               None
    minion     True  True         True          True     True       True         True                    False              1
    """

    data = StringIO(parameter_string_table)
    df = pd.read_csv(data, sep=r"\s+")
    run_properties = df.to_dict("records")

    return run_properties


@pytest.mark.parametrize("run_properties", build_run_properties())
def test_ont_transfer(create_dirs, run_properties, caplog):
    """Test the "taca analaysis ont-transfer" subcommand automation from
    start to finish for a variety of runs.
    """
    caplog.at_level(logging.INFO)

    # Create dir tree from fixture
    tmp = create_dirs

    # Create test config
    test_config_yaml = make_ONT_test_config(tmp)

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

    # Parametrized run
    create_ONT_run_dir(
        tmp,
        qc=run_properties["qc"],
        instrument=run_properties["instrument"],
        script_files=True,
        run_finished=run_properties["run_finished"],
        sync_finished=run_properties["sync_finished"],
        raw_dirs=run_properties["raw_dirs"],
        fastq_dirs=run_properties["fastq_dirs"],
        barcode_dirs=run_properties["barcode_dirs"],
        anglerfish_samplesheets=run_properties["anglerfish_samplesheets"],
    )

    # Start testing
    analysis_nanopore.ont_transfer(run_abspath=None, qc=False)
