#!/usr/bin/env python
import unittest
import logging
import filecmp
import mock
import os

from taca.analysis.analysis_nanopore import *
from taca.nanopore.minion import MinION
from taca.utils import config as conf


CONFIG = conf.load_yaml_config('data/taca_test_nanopore_cfg.yaml')

class TestNanoporeAnalysis(unittest.TestCase):
    def test_find_runs_to_process(self):
        """Find all expected nanopore runs to process."""
        expected_dirs = ["data/nanopore_data/run1/still_sequencing/20200101_1412_MN19414_AAU641_68125dc2",
                         "data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2",
                         "data/nanopore_data/run2/done_sequencing/20200102_1412_MN19414_AAU642_68125dc2",
                         "data/nanopore_data/run3/demultiplexing/20200103_1412_MN19414_AAU643_68125dc2",
                         "data/nanopore_data/run7/done_no_sample_sheet/20200107_1412_MN19417_AAU645_68125dc2",
                         "data/nanopore_data/run8/demux_failed/20200108_1412_MN19414_AAU648_68125dc2"]
        found_dirs = find_runs_to_process()
        self.assertEqual(sorted(found_dirs), sorted(expected_dirs))

    @mock.patch('taca.analysis.analysis_nanopore.os.path.isfile')
    @mock.patch('taca.nanopore.minion.MinION.start_nanoseq')
    def test_process_minion_run_start_analysis(self, mock_start, mock_isfile):
        """Start nanoseq analysis for minion."""
        nanoseq_sample_sheet = 'data/nanopore_data/run2/done_sequencing/20200102_1412_MN19414_AAU642_68125dc2/SQK-LSK109_sample_sheet.csv'
        anglerfish_sample_sheet = 'some/path'
        mock_isfile.return_value = True
        run_dir = 'data/nanopore_data/run2/done_sequencing/20200102_1412_MN19414_AAU642_68125dc2'
        minion_run = MinION(run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet)
        process_minion_run(minion_run)
        mock_start.assert_called_once()

    @mock.patch('taca.nanopore.minion.Nanopore.transfer_run')
    @mock.patch('taca.nanopore.minion.Nanopore.update_transfer_log')
    @mock.patch('taca.nanopore.minion.Nanopore.archive_run')
    @mock.patch('taca.analysis.analysis_nanopore.send_mail')
    def test_process_minion_run_transfer(self, mock_mail, mock_archive, mock_update, mock_transfer):
        """Start transfer of run directory."""
        mock_transfer.return_value = True
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        minion_run = MinION(run_dir, 'dummy/path', None)
        email_subject = ('Run successfully processed: 20200104_1412_MN19414_AAU644_68125dc2')
        email_message = 'Run 20200104_1412_MN19414_AAU644_68125dc2 has been analysed, transferred and archived successfully.'
        email_recipients = 'test@test.com'
        process_minion_run(minion_run)
        mock_mail.assert_called_once_with(email_subject, email_message, email_recipients)

    @mock.patch('taca.analysis.analysis_nanopore.send_mail')
    def test_process_minion_run_fail_analysis(self, mock_mail):
        """Send email to operator if nanoseq analysis failed."""
        run_dir = 'data/nanopore_data/run8/demux_failed/20200108_1412_MN19414_AAU648_68125dc2'
        minion_run = MinION(run_dir, None, None)
        process_minion_run(minion_run)
        email_subject = ('Analysis failed for run 20200108_1412_MN19414_AAU648_68125dc2')
        email_message = 'The nanoseq analysis failed for run {}.'.format(minion_run.run_id)
        email_recipients = 'test@test.com'
        mock_mail.assert_called_once_with(email_subject, email_message, email_recipients)
