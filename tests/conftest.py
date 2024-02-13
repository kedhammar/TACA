import os
import shutil
import tempfile

import pytest


@pytest.fixture
def create_dirs():
    """Create the bottom-level file-tree to be used for all tests:

        /var/folders/st/7720p4gn38vd042pfr7m2gh978ww82/T/tmp2viq7090
        ├── config
        │   └── Chromium_10X_indexes.txt
        ├── log
        │   ├── transfer_minion.tsv
        │   └── transfer_promethion.tsv
        ├── miarka
        │   ├── minion
        │   │   └── qc
        │   └── promethion
        ├── minknow_reports
        ├── ngi-nas-ns
        │   ├── NextSeq_data
        │   ├── NovaSeqXPlus_data
        │   ├── NovaSeq_data
        │   ├── minion_data
        │   ├── miseq_data
        │   ├── promethion_data
        │   └── samplesheets
        │       ├── NovaSeqXPlus
        │       └── anglerfish
        └── sequencing
            ├── MiSeq
            │   └── nosync
            ├── NextSeq
            │   └── nosync
            ├── NovaSeq
            │   └── nosync
            ├── NovaSeqXPlus
            │   └── nosync
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
    # Illumina
    os.makedirs(f"{tmp.name}/sequencing/MiSeq/nosync")
    os.makedirs(f"{tmp.name}/sequencing/NextSeq/nosync")
    os.makedirs(f"{tmp.name}/sequencing/NovaSeq/nosync")
    os.makedirs(f"{tmp.name}/sequencing/NovaSeqXPlus/nosync")
    # ONT
    os.makedirs(f"{tmp.name}/sequencing/promethion/nosync")
    os.makedirs(f"{tmp.name}/sequencing/minion/nosync")
    os.makedirs(f"{tmp.name}/sequencing/minion/qc/nosync")

    # Sequencing metadata
    # Illumina
    os.makedirs(f"{tmp.name}/ngi-nas-ns/miseq_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/NextSeq_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/NovaSeq_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/NovaSeqXPlus_data")
    # ONT
    os.makedirs(f"{tmp.name}/ngi-nas-ns/promethion_data")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/minion_data")

    # Samplesheets
    os.makedirs(f"{tmp.name}/ngi-nas-ns/samplesheets/anglerfish")
    os.makedirs(f"{tmp.name}/ngi-nas-ns/samplesheets/NovaSeqXPlus")

    # Misc. ONT dirs/files
    os.makedirs(f"{tmp.name}/minknow_reports")
    os.makedirs(f"{tmp.name}/log")
    open(f"{tmp.name}/log/transfer_promethion.tsv", "w").close()
    open(f"{tmp.name}/log/transfer_minion.tsv", "w").close

    # Analysis server destination dirs
    os.makedirs(f"{tmp.name}/miarka/promethion")
    os.makedirs(f"{tmp.name}/miarka/minion/qc")

    # Indexes
    os.makedirs(f"{tmp.name}/config")
    shutil.copy(
        "tests/data/Chromium_10X_indexes.txt",
        f"{tmp.name}/config/Chromium_10X_indexes.txt",
    )

    yield tmp

    tmp.cleanup()
