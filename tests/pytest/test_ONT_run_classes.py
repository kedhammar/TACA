import os
import tempfile

import pytest

from taca.nanopore import ONT_run_classes

# To check coverage, use
# pytest -s --cov=taca.nanopore.ONT_run_classes --cov-report term-missing -vv tests/pytest/test_ONT_run_classes.py


@pytest.fixture
def create_dirs():
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

    yield tmp

    tmp.cleanup()


def create_run_dir(
    tmp,
    instrument="promethion",
    instrument_position="1A",
    flowcell_id="TEST12345",
    data_dir=None,
):
    """Create a run directory according to specifications.

    ..
    └── {data_dir}
        └── 20240131_1702_{instrument_position}_{flowcell_id}_randomhash

    Return it's path.
    """
    if not data_dir:
        data_dir = f"{tmp.name}/sequencing/{instrument}_data"

    run_name = f"20240131_1702_{instrument_position}_{flowcell_id}_randomhash"
    run_path = f"{data_dir}/{run_name}"
    os.mkdir(run_path)

    return run_path


def test_ONT_user_run(create_dirs):
    tmp = create_dirs

    run_path = create_run_dir(tmp)
    run = ONT_run_classes.ONT_user_run(run_path)

    assert run.run_abspath == run_path
