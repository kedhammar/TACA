#!/usr/bin/env python
import unittest
import logging
import mock
import os

from taca.analysis.analysis_nanopore import *
from taca.nanopore.minion import MinIONqc
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
        nanopore_data_dir = CONFIG.get('nanopore_analysis').get('minion_qc_run').get('data_dir')
        skip_dirs = CONFIG.get('nanopore_analysis').get('minion_qc_run').get('ignore_dirs')
        found_dirs = find_minion_runs(nanopore_data_dir, skip_dirs)
        self.assertEqual(sorted(found_dirs), sorted(expected_dirs))

    @mock.patch('taca.analysis.analysis_nanopore.os.path.isfile')
    @mock.patch('taca.nanopore.minion.MinIONqc.start_nanoseq')
    def test_process_minion_run_start_analysis(self, mock_start, mock_isfile):
        """Start nanoseq analysis for minion."""
        nanoseq_sample_sheet = 'data/nanopore_data/run2/done_sequencing/20200102_1412_MN19414_AAU642_68125dc2/SQK-LSK109_sample_sheet.csv'
        anglerfish_sample_sheet = 'some/path'
        mock_isfile.return_value = True
        run_dir = 'data/nanopore_data/run2/done_sequencing/20200102_1412_MN19414_AAU642_68125dc2'
        minion_run = MinIONqc(run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet)
        process_minion_qc_run(minion_run)
        mock_start.assert_called_once()

    @mock.patch('taca.nanopore.minion.MinIONqc.copy_results_for_lims')
    @mock.patch('taca.nanopore.minion.Nanopore.transfer_run')
    @mock.patch('taca.nanopore.minion.Nanopore.update_transfer_log')
    @mock.patch('taca.nanopore.minion.Nanopore.archive_run')
    @mock.patch('taca.analysis.analysis_nanopore.send_mail')
    def test_process_minion_run_transfer(self, mock_mail, mock_archive, mock_update, mock_transfer, mock_cp):
        """Start transfer of run directory."""
        mock_transfer.return_value = True
        mock_cp.return_value = True
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        minion_run = MinIONqc(run_dir, 'dummy/path', None)
        email_subject = ('Run successfully processed: 20200104_1412_MN19414_AAU644_68125dc2')
        email_message = 'Run 20200104_1412_MN19414_AAU644_68125dc2 has been analysed, transferred and archived successfully.'
        email_recipients = 'test@test.com'
        process_minion_qc_run(minion_run)
        expected_calls = [mock.call('Anglerfish successfully processed run 20200104_1412_MN19414_AAU644_68125dc2', 
                           'Anglerfish has successfully finished for run 20200104_1412_MN19414_AAU644_68125dc2. Please finish the QC step in lims.', 
                           'test@test.com'),
                          mock.call('Run successfully processed: 20200104_1412_MN19414_AAU644_68125dc2', 
                           'Run 20200104_1412_MN19414_AAU644_68125dc2 has been analysed, transferred and archived successfully.', 
                           'test@test.com')]
        mock_mail.assert_has_calls(expected_calls)

    @mock.patch('taca.analysis.analysis_nanopore.send_mail')
    def test_process_minion_run_fail_analysis(self, mock_mail):
        """Send email to operator if nanoseq analysis failed."""
        run_dir = 'data/nanopore_data/run8/demux_failed/20200108_1412_MN19414_AAU648_68125dc2'
        minion_run = MinIONqc(run_dir, None, None)
        minion_run.qc_run = True
        process_minion_qc_run(minion_run)
        email_subject = ('Analysis failed for run 20200108_1412_MN19414_AAU648_68125dc2')
        email_message = 'The nanoseq analysis failed for run {}.'.format(minion_run.run_id)
        email_recipients = 'test@test.com'
        mock_mail.assert_called_once_with(email_subject, email_message, email_recipients)
