import json
import os
import re
import tempfile
from unittest.mock import Mock, call, mock_open, patch
from taca.nanopore import ONT_run_classes

import pytest

from taca.nanopore import instrument_transfer

# To check coverage, use
# pytest -s --cov=taca.nanopore.instrument_transfer --cov-report term-missing -vv test_instrument_transfer.py


@pytest.fixture
def create_dirs() -> tempfile.TemporaryDirectory:
    """Create the bottom-level file-tree to be used for all tests:

    tmp
    └── sequencing
        ├── promethion_data
        │   └── nosync
        └── minion_data
            ├── nosync
            └── qc
                └── nosync

    --> Return the the temporary directory object
    """
    tmp = tempfile.TemporaryDirectory()

    # Create the sequencing data dir structure
    os.makedirs(f"{tmp.name}/sequencing/promethion_data/nosync")
    os.makedirs(f"{tmp.name}/sequencing/minion_data/nosync")
    os.makedirs(f"{tmp.name}/sequencing/minion_data/qc/nosync")

    return tmp


def create_run_dir(
    tmp,
    instrument="promethion",
    instrument_position="1A",
    flowcell_id="TEST12345",
    data_dir=None,
):
    """Create a run directory.
    --> Return it's path.
    """
    if not data_dir:
        data_dir = f"{tmp.name}/sequencing/{instrument}_data/"

    run_name = f"20240131_1702_{instrument_position}_{flowcell_id}_randomhash"
    run_path = f"{data_dir}/{run_name}"
    os.mkdir(run_path)

    return run_path
