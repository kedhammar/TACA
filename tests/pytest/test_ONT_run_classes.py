import os
import tempfile
from unittest.mock import patch

import pytest
import yaml

from taca.nanopore import ONT_run_classes

# To check coverage, use
# pytest -s --cov=taca.nanopore.ONT_run_classes --cov-report term-missing -vv tests/pytest/test_ONT_run_classes.py


def make_test_config(tmp):
    test_config_yaml_string = f"""nanopore_analysis:
    run_types:
        user_run:
        data_dirs:
            - {tmp.name}/sequencing/promethion
            - {tmp.name}/sequencing/minion
        ignore_dirs:
            - 'nosync'
            - 'qc'
        instruments:
            promethion:
                transfer_log: {tmp.name}/log/transfer_promethion.tsv
                archive_dir: {tmp.name}/sequencing/promethion/nosync
                metadata_dir: {tmp.name}/ngi-nas-ns/promethion_data
                destination: {tmp.name}/miarka/promethion/
            minion:
                transfer_log: /{tmp.name}/log/transfer_minion.tsv
                archive_dir: {tmp.name}/sequencing/minion/nosync
                metadata_dir: {tmp.name}/ngi-nas-ns/minion_data
                destination: {tmp.name}/miarka/minion/
        qc_run:
            data_dirs:
                - {tmp.name}/ngi_data/sequencing/minion/qc
        ignore_dirs:
            - 'nosync'
        instruments:
            minion:
                transfer_log: {tmp.name}/log/transfer_minion_qc.tsv
                archive_dir: {tmp.name}/sequencing/minion/qc/nosync
                metadata_dir: {tmp.name}/ngi-nas-ns/minion_data/qc
                destination: {tmp.name}/miarka/minion/qc
        anglerfish:
            anglerfish_samplesheets_dir: /srv/ngi-nas-ns/samplesheets/anglerfish
            anglerfish_path: ~/miniconda3/envs/anglerfish/bin/anglerfish
    minknow_reports_dir: {tmp.name}/minknow_reports/
    analysis_server:
        host: miarka1.uppmax.uu.se
        port:
        user: funk_903
    rsync_options:
        '-Lav': None
        '--chown': ':ngi2016003'
        '--chmod': 'Dg+s,g+rw'
        '-r': None
        '--exclude': ['work']"""

    test_config_yaml = yaml.safe_load(test_config_yaml_string)

    return test_config_yaml


@pytest.fixture
def create_dirs():
    """Create the bottom-level file-tree to be used for all tests:

    tmp
    ├── minknow_reports
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

    # CREATE DIR STRUCTURE

    # Sequencing data
    os.makedirs(f"{tmp.name}/sequencing/promethion_data/nosync")
    os.makedirs(f"{tmp.name}/sequencing/minion_data/nosync")
    os.makedirs(f"{tmp.name}/sequencing/minion_data/qc/nosync")

    # Non-sensitive metadata
    os.makedirs(f"{tmp.name}/ngi-nas-ns/promethion_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/minion_data")

    # Reports for GenStat
    os.makedirs(f"{tmp.name}/minknow_reports")

    # Logs
    os.makedirs(f"{tmp.name}/log")
    open(f"{tmp.name}/log/transfer_promethion.tsv", "w").close()
    open(f"{tmp.name}/log/transfer_minion.tsv", "w").close

    # Analysis server destination dirs
    os.makedirs(f"{tmp.name}/miarka/promethion")
    os.makedirs(f"{tmp.name}/miarka/minion/qc")

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

    # Mock the CONFIG
    mock_config = patch("taca.utils.config.CONFIG", new=make_test_config(tmp))
    mock_config.start()

    run_path = create_run_dir(tmp)
    run = ONT_run_classes.ONT_user_run(run_path)

    assert run.run_abspath == run_path
