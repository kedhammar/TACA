#!/usr/bin/env python

import os
import shutil
import tempfile
import unittest
import csv
import json
import mock
from datetime import datetime

from taca.analysis.analysis import *
from taca.illumina.Runs import Run
from taca.illumina.HiSeqX_Runs import HiSeqX_Run
from taca.utils import config as conf


# This is only run if TACA is called from the CLI, as this is a test, we need to
# call it explicitely
CONFIG = conf.load_yaml_config('data/taca_test_cfg.yaml')

def processing_status(run_dir):
    demux_dir = os.path.join(run_dir, 'Demultiplexing')
    if not os.path.exists(demux_dir):
        return 'TO_START'
    elif os.path.exists(os.path.join(demux_dir, 'Stats', 'Stats.json')):
        return 'COMPLETED'
    else:
        return 'IN_PROGRESS'

class TestTracker(unittest.TestCase):
    """ analysis.py script tests
    """
    @classmethod
    def setUpClass(self):
        """ Creates the following directory tree for testing purposes:

        tmp/
        |__ 141124_ST-COMPLETED_01_AFCIDXX
        |   |__ RunInfo.xml
        |   |__ Demultiplexing
        |   |   |__ Stats
        |   |       |__ DemultiplexingStats.xml
        |   |__ RTAComplete.txt
        |__ 141124_ST-INPROGRESS_02_AFCIDXX
        |   |__ RunInfo.xml
        |   |__ Demultiplexing
        |   |__ Demultiplexing_0
        |   |__ Demultiplexing_1
        |   |__ SampleSheet_0.csv
        |   |__ SampleSheet_1.csv
        |   |__ RTAComplete.txt
        |__ 141124_ST-INPROGRESSDONE_02_AFCIDXX
        |   |__ RunInfo.xml
        |   |__ Demultiplexing
        |   |__ Demultiplexing_0
        |   |   |__Stats
        |   |      |__ DemultiplexingStats.xml
        |   |__ Demultiplexing_1
        |   |   |__Stats
        |   |      |__ DemultiplexingStats.xml
        |   |__ SampleSheet_0.csv
        |   |__ SampleSheet_1.csv
        |   |__ RTAComplete.txt
        |__ 141124_ST-RUNNING_03_AFCIDXX
        |   |__ RunInfo.xml
        |__ 141124_ST-TOSTART_04_FCIDXXX
        |   |__ RunInfo.xml
        |   |__ RTAComplete.txt
        |__
        """
        self.tmp_dir = os.path.join(tempfile.mkdtemp(), 'tmp')
        self.transfer_file = os.path.join(self.tmp_dir, 'transfer.tsv')

        running = os.path.join(self.tmp_dir, '141124_ST-RUNNING1_03_AFCIDXX')
        to_start = os.path.join(self.tmp_dir, '141124_ST-TOSTART1_04_FCIDXXX')
        in_progress = os.path.join(self.tmp_dir, '141124_ST-INPROGRESS1_02_AFCIDXX')
        in_progress_done = os.path.join(self.tmp_dir, '141124_ST-INPROGRESSDONE1_02_AFCIDXX')
        completed = os.path.join(self.tmp_dir, '141124_ST-COMPLETED1_01_AFCIDXX')
        finished_runs = [to_start, in_progress, in_progress_done, completed]

        # Create runs directory structure
        os.makedirs(self.tmp_dir)
        os.makedirs(running)
        os.makedirs(to_start)
        os.makedirs(os.path.join(in_progress, 'Demultiplexing'))
        os.makedirs(os.path.join(in_progress, 'Demultiplexing_0'))
        os.makedirs(os.path.join(in_progress, 'Demultiplexing_1'))
        os.makedirs(os.path.join(in_progress_done, 'Demultiplexing'))
        os.makedirs(os.path.join(in_progress_done, 'Demultiplexing_0/Stats'))
        os.makedirs(os.path.join(in_progress_done, 'Demultiplexing_1/Stats'))
        os.makedirs(os.path.join(completed, 'Demultiplexing', 'Stats'))

        # Create files indicating that the run is finished
        for run in finished_runs:
            open(os.path.join(run, 'RTAComplete.txt'), 'w').close()

        # Create sample sheets for running demultiplexing
        open(os.path.join(in_progress, 'SampleSheet_0.csv'), 'w').close()
        open(os.path.join(in_progress, 'SampleSheet_1.csv'), 'w').close()
        open(os.path.join(in_progress_done, 'SampleSheet_0.csv'), 'w').close()
        open(os.path.join(in_progress_done, 'SampleSheet_1.csv'), 'w').close()

        # Create files indicating that demultiplexing is ongoing
        open(os.path.join(in_progress_done, 'Demultiplexing_0', 'Stats', 'DemultiplexingStats.xml'), 'w').close()
        open(os.path.join(in_progress_done, 'Demultiplexing_1', 'Stats', 'DemultiplexingStats.xml'), 'w').close()

        # Create files indicating that the preprocessing is done
        open(os.path.join(completed, 'Demultiplexing', 'Stats', 'DemultiplexingStats.xml'), 'w').close()
        with open(os.path.join(completed, 'Demultiplexing', 'Stats', 'Stats.json'), 'w') as stats_json:
            json.dump({"silly": 1}, stats_json)

        # Create transfer file and add the completed run
        with open(self.transfer_file, 'w') as f:
            tsv_writer = csv.writer(f, delimiter='\t')
            tsv_writer.writerow([os.path.basename(completed), str(datetime.now())])

        # Move sample RunInfo.xml file to every run directory
        for run in [running, to_start, in_progress, in_progress_done, completed]:
            shutil.copy('data/RunInfo.xml', run)
            shutil.copy('data/runParameters.xml', run)

        # Create run objects
        self.running = HiSeqX_Run(os.path.join(self.tmp_dir,
                                               '141124_ST-RUNNING1_03_AFCIDXX'),
                                  CONFIG["analysis"]["HiSeqX"])
        self.to_start = Run(os.path.join(self.tmp_dir,
                                         '141124_ST-TOSTART1_04_FCIDXXX'),
                            CONFIG["analysis"]["HiSeqX"])
        self.in_progress = HiSeqX_Run(os.path.join(self.tmp_dir,
                                            '141124_ST-INPROGRESS1_02_AFCIDXX'),
                               CONFIG["analysis"]["HiSeqX"])
        self.in_progress_done = HiSeqX_Run(os.path.join(self.tmp_dir,
                                            '141124_ST-INPROGRESSDONE1_02_AFCIDXX'),
                               CONFIG["analysis"]["HiSeqX"])
        self.completed = Run(os.path.join(self.tmp_dir,
                                          '141124_ST-COMPLETED1_01_AFCIDXX'),
                             CONFIG["analysis"]["HiSeqX"])
        self.finished_runs = [self.to_start, self.in_progress, self.completed]
        self.transfer_file = os.path.join(self.tmp_dir, 'transfer.tsv')

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.tmp_dir)

    def test_is_sequencing_done(self):
        """ Is finished should be True only if "RTAComplete.txt" file is present...
        """
        self.assertFalse(self.running._is_sequencing_done())
        self.assertTrue(all(map(lambda run: run._is_sequencing_done, self.finished_runs)))

    def test_get_run_status(self):
        """ Status of the processing depends on the generated files
        """
        self.assertEqual('SEQUENCING', self.running.get_run_status())
        self.assertEqual('TO_START', self.to_start.get_run_status())
        self.assertEqual('IN_PROGRESS', self.in_progress.get_run_status())
        self.assertEqual('COMPLETED', self.completed.get_run_status())

    def test_is_transferred(self):
        """ is_transferred should rely on the info in transfer.tsv
        """
        self.assertTrue(self.completed.is_transferred(self.transfer_file))
        self.assertFalse(self.running.is_transferred(self.transfer_file))
        self.assertFalse(self.to_start.is_transferred(self.transfer_file))
        self.assertFalse(self.in_progress.is_transferred( self.transfer_file))

    @mock.patch('taca.illumina.HiSeqX_Runs.HiSeqX_Run._aggregate_demux_results')
    def test_check_run_status_done(self, mock_aggregate_demux_results):
        """ check_run_status should recognize if a demultiplexing run is finished or not
        """
        # _aggregate_demux_results is only run if the demultiplexing is done
        self.in_progress.check_run_status()
        mock_aggregate_demux_results.assert_not_called()
        self.in_progress_done.check_run_status()
        mock_aggregate_demux_results.assert_called_once()

    def test_get_run_type(self):
        """ Return runtype if set"""
        self.assertEqual('NGI-RUN', self.running.get_run_type())

    def test_get_demux_folder(self):
        """ Return name of demux folder if set"""
        self.assertEqual('Demultiplexing', self.running._get_demux_folder())

    def test_get_samplesheet(self):
        """ Return location of sample sheet"""
        self.assertEqual('data/2014/FCIDXX.csv', self.running._get_samplesheet())

    def test_is_demultiplexing_done(self):
        """ Return true if Stats.json exists, else false"""
        self.assertFalse(self.in_progress._is_demultiplexing_done())
        self.assertTrue(self.completed._is_demultiplexing_done())

    def test_is_demultiplexing_started(self):
        """ Return true if demux folder exists, else false"""
        self.assertTrue(self.in_progress._is_demultiplexing_started())
        self.assertFalse(self.to_start._is_demultiplexing_started())


