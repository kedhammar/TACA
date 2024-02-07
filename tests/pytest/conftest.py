import os
import tempfile

import pytest


@pytest.fixture
def create_dirs():
    """Create the bottom-level file-tree to be used for all tests:

    tmp
    ├── log
    │   ├── transfer_minion.tsv
    │   └── transfer_promethion.tsv
    ├── miarka
    │   ├── minion
    │   │   └── qc
    │   └── promethion
    ├── minknow_reports
    ├── ngi-nas-ns
    │   ├── minion_data
    │   ├── promethion_data
    │   └── samplesheets
    │       └── anglerfish
    └── sequencing
        ├── minion
        │   ├── nosync
        │   └── qc
        │       └── nosync
        └── promethion
            └── nosync

    --> Return the the temporary directory object
    """
    tmp = tempfile.TemporaryDirectory()

    # CREATE DIR STRUCTURE

    # Sequencing data
    os.makedirs(f"{tmp.name}/sequencing/promethion/nosync")
    os.makedirs(f"{tmp.name}/sequencing/minion/nosync")
    os.makedirs(f"{tmp.name}/sequencing/minion/qc/nosync")

    # Non-sensitive metadata
    os.makedirs(f"{tmp.name}/ngi-nas-ns/promethion_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/minion_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/samplesheets/anglerfish")

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
