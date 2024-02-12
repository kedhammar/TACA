import importlib
import subprocess
from unittest.mock import patch

import yaml

from taca.analysis import analysis


def make_illumina_test_config(tmp):
    test_config_yaml_string = f"""mail: 
    recipients: mock
statusdb: mock
analysis:
    MiSeq:
        QC:
            max_number_undetermined_reads_simple_lane: 5000000
            max_percentage_undetermined_indexes_simple_lane: 5
            max_percentage_undetermined_indexes_NoIndex_lane: 30
            max_frequency_most_represented_und_index_NoIndex_lane: 40
        bclconvert:
            Version: bcl-convert Version 00.000.000.4.2.4
            bin: mock
            tenX_index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            smartseq_index_path: {tmp.name}/config/Smart-seq3_v1.5.csv
            index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            options:
                common:
                    - bcl-sampleproject-subdirectories: true
                    - fastq-gzip-compression-level: 4
                    - output-directory: Demultiplexing
                    - output-legacy-stats: true
                    - sample-name-column-enabled: true
                    - shared-thread-odirect-output: true
                    - strict-mode: false
                    - force
                SMARTSEQ:
                    - no-lane-splitting: true
            settings:
                common:
                    - MinimumTrimmedReadLength: 0
                    - MaskShortReads: 0
                10X_SINGLE:
                    - CreateFastqForIndexReads: 1
                IDT_UMI:
                    - CreateFastqForIndexReads: 1
                    - TrimUMI: 0
                SMARTSEQ:
                    - CreateFastqForIndexReads: 1
                NOINDEX:
                    - CreateFastqForIndexReads: 1
                    - BarcodeMismatchesIndex1: 0
                    - BarcodeMismatchesIndex2: 0
                short_single_index:
                    - BarcodeMismatchesIndex1: 0
                    - BarcodeMismatchesIndex2: 0
        bcl2fastq:
            Version: bcl2fastq_v2.20.0.422
            bin: mock
            index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            tenX_index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            smartseq_index_path: {tmp.name}/config/Smart-seq3_v1.5.csv
            options:
                common:
                    - output-dir: Demultiplexing
                    - loading-threads: 1
                    - processing-threads: 8
                    - writing-threads: 1
                    - minimum-trimmed-read-length: 0
                    - mask-short-adapter-reads: 0
                    - ignore-missing-positions
                    - ignore-missing-controls
                    - ignore-missing-filter
                    - ignore-missing-bcls
                10X_SINGLE:
                    - create-fastq-for-index-reads
                SMARTSEQ:
                    - create-fastq-for-index-reads
                    - no-lane-splitting
                NOINDEX:
                    - create-fastq-for-index-reads
                    - barcode-mismatches: 0
                short_single_index:
                    - barcode-mismatches: 0
        samplesheets_dir: ''
        analysis_server:
            host: miarka1.uppmax.uu.se
            port:
            user: funk_903
            sync:
                data_archive: ""
                include:
                    - "*.xml"
                    - "*.htm"
                    - "*.html"
                    - "*.csv"
                    - "*.fastq.gz"
    NovaSeq:
        bclconvert:
            Version: bcl-convert Version 00.000.000.4.2.4
            bin: mock
            tenX_index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            smartseq_index_path: {tmp.name}/config/Smart-seq3_v1.5.csv
            index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            options:
                common:
                    - bcl-sampleproject-subdirectories: true
                    - fastq-gzip-compression-level: 4
                    - output-directory: Demultiplexing
                    - output-legacy-stats: true
                    - sample-name-column-enabled: true
                    - shared-thread-odirect-output: true
                    - strict-mode: false
                    - force
                SMARTSEQ:
                    - no-lane-splitting: true
            settings:
                common:
                    - MinimumTrimmedReadLength: 0
                    - MaskShortReads: 0
                10X_SINGLE:
                    - CreateFastqForIndexReads: 1
                IDT_UMI:
                    - CreateFastqForIndexReads: 1
                    - TrimUMI: 0
                SMARTSEQ:
                    - CreateFastqForIndexReads: 1
                NOINDEX:
                    - CreateFastqForIndexReads: 1
                    - BarcodeMismatchesIndex1: 0
                    - BarcodeMismatchesIndex2: 0
                short_single_index:
                    - BarcodeMismatchesIndex1: 0
                    - BarcodeMismatchesIndex2: 0
        bcl2fastq:
            Version: bcl2fastq_v2.20.0.422
            bin: mock
            index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            tenX_index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            smartseq_index_path: {tmp.name}/config/Smart-seq3_v1.5.csv
            options:
                common:
                    - output-dir: Demultiplexing
                    - loading-threads: 1
                    - processing-threads: 8
                    - writing-threads: 1
                    - minimum-trimmed-read-length: 0
                    - mask-short-adapter-reads: 0
                    - ignore-missing-positions
                    - ignore-missing-controls
                    - ignore-missing-filter
                    - ignore-missing-bcls
                10X_SINGLE:
                    - create-fastq-for-index-reads
                SMARTSEQ:
                    - create-fastq-for-index-reads
                    - no-lane-splitting
                NOINDEX:
                    - create-fastq-for-index-reads
                    - barcode-mismatches: 0
                short_single_index:
                    - barcode-mismatches: 0
        samplesheets_dir: {tmp.name}/ngi-nas-ns/samplesheets/novaseq
        analysis_server:
            host: miarka1.uppmax.uu.se
            port:
            user: funk_903
            sync:
                data_archive: ""
                include:
                    - "*.xml"
                    - "*.htm"
                    - "*.html"
                    - "*.fastq.gz"
                    - "*.csv"
    NovaSeqXPlus:
        bclconvert:
            Version: bcl-convert Version 00.000.000.4.2.4
            bin: mock
            tenX_index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            smartseq_index_path: {tmp.name}/config/Smart-seq3_v1.5.csv
            index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            options:
                common:
                    - bcl-sampleproject-subdirectories: true
                    - fastq-gzip-compression-level: 4
                    - output-directory: Demultiplexing
                    - output-legacy-stats: true
                    - sample-name-column-enabled: true
                    - shared-thread-odirect-output: true
                    - strict-mode: false
                    - force
                SMARTSEQ:
                    - no-lane-splitting: true
            settings:
                common:
                    - MinimumTrimmedReadLength: 0
                    - MaskShortReads: 0
                10X_SINGLE:
                    - CreateFastqForIndexReads: 1
                IDT_UMI:
                    - CreateFastqForIndexReads: 1
                    - TrimUMI: 0
                SMARTSEQ:
                    - CreateFastqForIndexReads: 1
                NOINDEX:
                    - CreateFastqForIndexReads: 1
                    - BarcodeMismatchesIndex1: 0
                    - BarcodeMismatchesIndex2: 0
                short_single_index:
                    - BarcodeMismatchesIndex1: 0
                    - BarcodeMismatchesIndex2: 0
        bcl2fastq:
            Version: bcl2fastq_v2.20.0.422
            bin: mock
            tenX_index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            smartseq_index_path: {tmp.name}/config/Smart-seq3_v1.5.csv
            index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            options:
                common:
                    - output-dir: Demultiplexing
                    - loading-threads: 1
                    - processing-threads: 8
                    - writing-threads: 1
                    - minimum-trimmed-read-length: 0
                    - mask-short-adapter-reads: 0
                    - ignore-missing-positions
                    - ignore-missing-controls
                    - ignore-missing-filter
                    - ignore-missing-bcls
                10X_SINGLE:
                    - create-fastq-for-index-reads
                SMARTSEQ:
                    - create-fastq-for-index-reads
                    - no-lane-splitting
                NOINDEX:
                    - create-fastq-for-index-reads
                    - barcode-mismatches: 0
                short_single_index:
                    - barcode-mismatches: 0
        samplesheets_dir: {tmp.name}/ngi-nas-ns/samplesheets/NovaSeqXPlus
        analysis_server:
            host: miarka1.uppmax.uu.se
            port:
            user: funk_903
            sync:
                data_archive: ""
                include:
                    - "*.xml"
                    - "*.htm"
                    - "*.html"
                    - "*.fastq.gz"
                    - "*.csv"
    NextSeq:
        bclconvert:
            Version: bcl-convert Version 00.000.000.4.2.4
            bin: mock
            tenX_index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            smartseq_index_path: {tmp.name}/config/Smart-seq3_v1.5.csv
            index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            options:
                common:
                    - bcl-sampleproject-subdirectories: true
                    - fastq-gzip-compression-level: 4
                    - output-directory: Demultiplexing
                    - output-legacy-stats: true
                    - sample-name-column-enabled: true
                    - shared-thread-odirect-output: true
                    - strict-mode: false
                    - force
                SMARTSEQ:
                    - no-lane-splitting: true
            settings:
                common:
                    - MinimumTrimmedReadLength: 0
                    - MaskShortReads: 0
                10X_SINGLE:
                    - CreateFastqForIndexReads: 1
                IDT_UMI:
                    - CreateFastqForIndexReads: 1
                    - TrimUMI: 0
                SMARTSEQ:
                    - CreateFastqForIndexReads: 1
                NOINDEX:
                    - CreateFastqForIndexReads: 1
                    - BarcodeMismatchesIndex1: 0
                    - BarcodeMismatchesIndex2: 0
                short_single_index:
                    - BarcodeMismatchesIndex1: 0
                    - BarcodeMismatchesIndex2: 0
        bcl2fastq:
            Version: bcl2fastq_v2.20.0.422
            bin: {tmp.name}/src/bcl2fastq_v2.20.0.422/bin/bcl2fastq
            index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            tenX_index_path: {tmp.name}/config/Chromium_10X_indexes.txt
            smartseq_index_path: {tmp.name}/config/Smart-seq3.csv
            options:
                common:
                    - output-dir: Demultiplexing
                    - loading-threads: 1
                    - processing-threads: 8
                    - writing-threads: 1
                    - minimum-trimmed-read-length: 0
                    - mask-short-adapter-reads: 0
                    - ignore-missing-positions
                    - ignore-missing-controls
                    - ignore-missing-filter
                    - ignore-missing-bcls
                10X_SINGLE:
                    - create-fastq-for-index-reads
                SMARTSEQ:
                    - create-fastq-for-index-reads
                    - no-lane-splitting
                NOINDEX:
                    - create-fastq-for-index-reads
                    - barcode-mismatches: 0
                short_single_index:
                    - barcode-mismatches: 0
        samplesheets_dir: {tmp.name}/ngi-nas-ns/samplesheets/nextseq
        analysis_server:
            host: miarka1.uppmax.uu.se
            port:
            user: funk_903
            sync:
                data_archive: ""
                include:
                    - "*.xml"
                    - "*.htm"
                    - "*.html"
                    - "*.fastq.gz"
                    - "*.csv"

   # Directories where TACA will look for new data to be processed
    data_dirs:
        - {tmp.name}/ngi_data/sequencing/MiSeq
        - {tmp.name}/ngi_data/sequencing/NovaSeq
        - {tmp.name}/ngi_data/sequencing/NovaSeqXPlus
        - {tmp.name}/ngi_data/sequencing/NextSeq/Runs
    # finished_dir is the direcory suffix that TACA will use to store flowcells after being
    # processed and transferred. The end finished_dir directories will be
    # data_dir/finished_dir for data_dir in data_dirs
    finished_dir: nosync
    # status_dir indicates the directory where TACA will store the CSV file used
    # to keep track of which flowcells have already been processed and transferred
    status_dir: {tmp.name}/log/
    mfs_path:
        miseq: {tmp.name}/ngi-nas-ns/miseq_data
        novaseq: {tmp.name}/ngi-nas-ns/NovaSeq_data
        novaseqxplus: {tmp.name}/ngi-nas-ns/NovaSeqXPlus_data
        nextseq: {tmp.name}/ngi-nas-ns/NextSeq_data
    deliver_runfolder:
        analysis_server:
            host: miarka1.uppmax.uu.se
            port:
            user: funk_903
        destination: /"""

    test_config_yaml = yaml.safe_load(test_config_yaml_string)

    return test_config_yaml


def create_illumina_run_dir(
    tmp,
):
    """Create a run directory according to specifications.

    Return it's path.
    """


def test_analysis(create_dirs):
    tmp = create_dirs

    test_config_yaml = make_illumina_test_config(tmp)

    assert test_config_yaml
