#!/usr/bin/env python

#import os
#import shutil
#import tempfile
import unittest
#import csv
#import json

#from datetime import datetime

from taca.analysis.analysis_nanopore import *
from taca.utils import config as conf


# This is only run if TACA is called from the CLI, as this is a test, we need to
# call it explicitely
CONFIG = conf.load_yaml_config('data/taca_test_nanopore_cfg.yaml')

class TestNanoporeAnalysis(unittest.TestCase):

    def test_find_runs_to_process(self):
        ''' Find all expected nanopore runs to process
        '''
        expected_dirs = ["data/nanopore_data/run1/run1/run1_sequencing_rna",
                         "data/nanopore_data/run4/run4/run4_done_demuxing_rna",
                         "data/nanopore_data/run2/run2/run2_done_sequencing_rna",
                         "data/nanopore_data/run3/run3/run3_demuxing_rna",
                         "data/nanopore_data/run7/run7/run7_done_no_sample_sheet",
                         "data/nanopore_data/run8/run8/run8_demux_failed_rna"]
        found_dirs = find_runs_to_process()
        self.assertItemsEqual(found_dirs, expected_dirs) # Prob not going to work since found dirs are returned in random order...

    def test_check_exit_status(self):
        ''' Check nanoseq exit status from file
        '''
        self.assertTrue(check_exit_status("data/nanopore_data/run4/run4/run4_done_demuxing_rna/.exitcode_for_nanoseq"))
        self.assertFalse(check_exit_status("data/nanopore_data/run8/run8/run8_demux_failed_rna/.exitcode_for_nanoseq"))

    def test_is_not_transferred(self):
        ''' Check if nanopore run has been transferred
        '''
        self.assertTrue(is_not_transferred("run4_done_demuxing_rna", "data/nanopore_data/transfer.tsv"))
        self.assertFalse(is_not_transferred("run5_completed_and_moved", "data/nanopore_data/transfer.tsv"))
