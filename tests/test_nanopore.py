#!/usr/bin/env python
import unittest
import mock
import filecmp
import os
import subprocess

from taca.nanopore.nanopore import Nanopore
from taca.nanopore.promethion import PromethION
from taca.nanopore.minion import MinION
from taca.utils import config as conf

CONFIG = conf.load_yaml_config('data/taca_test_nanopore_cfg.yaml')

class TestNanopore(unittest.TestCase):
    """Test Nanopore class"""
    def test_is_not_transferred(self):
        """Check if nanopore run has been transferred."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        np_run = Nanopore(run_dir)
        self.assertTrue(np_run.is_not_transferred())
        run_dir_transf = 'data/nanopore_data/run4/done_demuxing/20200105_1412_MN19414_AAU645_68125dc2'
        np_run_transf = Nanopore(run_dir_transf)
        self.assertFalse(np_run_transf.is_not_transferred())

    @mock.patch('taca.nanopore.nanopore.RsyncAgent')
    def test_transfer_run(self, mock_rsync):
        """Start rsync of finished run."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        np_run = Nanopore(run_dir)
        np_run.transfer_run()
        rsync_opts = {'-Lav': None,
                      '--chown': ':ngi2016003',
                      '--chmod' : 'Dg+s,g+rw',
                      '-r' : None,
                      '--exclude' : 'work'}
        mock_rsync.assert_called_with(run_dir,
                                      dest_path='some_dir',
                                      remote_host='some_host',
                                      remote_user='some_user',
                                      validate=False,
                                      opts=rsync_opts)

    @mock.patch('taca.nanopore.nanopore.shutil.move')
    def test_archive_run(self, mock_move):
        """Move directory to archive."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        np_run = Nanopore(run_dir)
        np_run.archive_run()
        mock_move.assert_called_once()

    
class TestMinION(unittest.TestCase):
    """Test MinION class"""

    def test_get_original_samplesheet(self):
        """Get location of lims sample sheet."""
        run_dir = 'data/nanopore_data/run2/done_sequencing/20200102_1412_MN19414_AAU642_68125dc2'
        run = MinION(run_dir, None, None)
        run._get_original_samplesheet()
        expected_sample_sheet = 'data/nanopore_samplesheets/2020/QC_SQK-LSK109_AAU642_Samplesheet_22-594126.csv'
        self.assertEqual(run.lims_samplesheet, expected_sample_sheet)

    def test_parse_samplesheet(self):
        """Make nanoseq sample sheet from lims sample sheet."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        run = MinION(run_dir, None, None)
        run.lims_samplesheet = 'data/nanopore_samplesheets/2020/DELIVERY_SQK-LSK109_AAU644_Samplesheet_24-594126.csv'
        run._parse_samplesheet()
        self.assertTrue(filecmp.cmp(run.nanoseq_sample_sheet, 'data/nanopore_samplesheets/expected/SQK-LSK109_sample_sheet.csv'))
        self.assertTrue(filecmp.cmp(run.anglerfish_sample_sheet, 'data/nanopore_samplesheets/expected/anglerfish_sample_sheet.csv'))

    @mock.patch('taca.nanopore.minion.MinION._get_flowcell_product_code')
    @mock.patch('taca.nanopore.minion.MinION._is_multiplexed')
    @mock.patch('taca.nanopore.minion.subprocess.Popen')
    def test_start_analysis_pipeline_multiplexed(self, mock_popen, mock_is_multiplexed, mock_get_fc_code):
        """Submit detached nanoseq job for multiplexed data."""
        mock_get_fc_code.return_value = 'FLO-FLG001'
        mock_is_multiplexed.return_value = True
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        sample_sheet = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2/SQK-LSK109_sample_sheet.csv'
        run = MinION(run_dir, sample_sheet, None)
        run.start_nanoseq()
        expected_parameters = ('nextflow run nf-core/nanoseq --input ' + sample_sheet
                               + ' --input_path ' + os.path.join(run_dir, 'fast5')
                               + ' --outdir ' + os.path.join(run_dir, 'nanoseq_output')
                               + ' --flowcell FLO-FLG001'
                               + ' --guppy_gpu'
                               + ' --skip_alignment'
                               + ' --kit SQK-LSK109'
                               + ' --max_cpus 6'
                               + ' --max_memory 20.GB'
                               + ' --barcode_kit EXP-NBD104'
                               + ' -profile singularity; echo $? > .exitcode_for_nanoseq')
        mock_popen.assert_called_once_with(expected_parameters, stdout=subprocess.PIPE, shell=True, cwd=run_dir)

    @mock.patch('taca.nanopore.minion.MinION._get_flowcell_product_code')
    @mock.patch('taca.nanopore.minion.MinION._is_multiplexed')
    @mock.patch('taca.nanopore.minion.subprocess.Popen')
    def test_start_analysis_pipeline_not_multiplexed(self, mock_popen, mock_is_multiplexed, mock_get_fc_code):
        """Submit detached nanoseq job for non multiplexed data."""
        mock_get_fc_code.return_value = 'FLO-FLG001'
        mock_is_multiplexed.return_value = False
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        sample_sheet = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2/SQK-LSK109_sample_sheet.csv'
        run = MinION(run_dir, sample_sheet, None)
        run.start_nanoseq()
        expected_parameters = ('nextflow run nf-core/nanoseq --input ' + sample_sheet
                               + ' --input_path ' + os.path.join(run_dir, 'fast5')
                               + ' --outdir ' + os.path.join(run_dir, 'nanoseq_output')
                               + ' --flowcell FLO-FLG001'
                               + ' --guppy_gpu'
                               + ' --skip_alignment'
                               + ' --kit SQK-LSK109'
                               + ' --max_cpus 6'
                               + ' --max_memory 20.GB'
                               + ' -profile singularity; echo $? > .exitcode_for_nanoseq')
        mock_popen.assert_called_once_with(expected_parameters, stdout=subprocess.PIPE, shell=True, cwd=run_dir)

    def test_get_flowcell_product_code(self):
        """Get flowcell product code from report.md."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        run = MinION(run_dir, None, None)
        got_id = run._get_flowcell_product_code()
        expected_id = 'FLO-FLG001'
        self.assertEqual(got_id, expected_id)

    def test_is_multiplexed(self):
        """Return True if run is multiplexed, else False."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        multiplexed_sample_sheet = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2/SQK-LSK109_sample_sheet.csv'
        non_multiplexed_sample_sheet = 'data/nanopore_data/run3/demultiplexing/20200103_1412_MN19414_AAU643_68125dc2/SQK-LSK109_AAU643_sample_sheet.csv'
        multiplexed_run = MinION(run_dir, multiplexed_sample_sheet, None)
        non_multiplexed_run = MinION(run_dir, non_multiplexed_sample_sheet, None)
        self.assertTrue(multiplexed_run._is_multiplexed())
        self.assertFalse(non_multiplexed_run._is_multiplexed())

    def test_get_barcode_kit(self):
        """Return EXP-NBD104 or EXP-NBD114 barcode kit based on sample sheet."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        sample_sheet_104 = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2/SQK-LSK109_sample_sheet.csv'
        run_104 = MinION(run_dir, sample_sheet_104, None)
        got_kit_104 = run_104._get_barcode_kit()
        
        sample_sheet_114 = 'data/nanopore_data/run8/demux_failed/20200108_1412_MN19414_AAU648_68125dc2/SQK-LSK109_sample_sheet.csv'
        run_114 = MinION(run_dir, sample_sheet_114, None)
        got_kit_114 = run_114._get_barcode_kit()
        self.assertEqual(got_kit_104, 'EXP-NBD104')
        self.assertEqual(got_kit_114, 'EXP-NBD114')

    def test_check_exit_status(self):
        """Check nanoseq exit status from file."""
        run_dir_success = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        success_run = MinION(run_dir_success, None, None)
        self.assertTrue(success_run.check_exit_status('data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2/.exitcode_for_nanoseq'))
        run_dir_fail = 'data/nanopore_data/run8/demux_failed/20200108_1412_MN19414_AAU648_68125dc2'
        fail_run = MinION(run_dir_fail, None, None)
        self.assertFalse(fail_run.check_exit_status('data/nanopore_data/run8/demux_failed/20200108_1412_MN19414_AAU648_68125dc2/.exitcode_for_nanoseq'))

    @mock.patch('taca.nanopore.minion.os.makedirs')
    @mock.patch('taca.nanopore.minion.subprocess.Popen')
    def test_start_anglerfish(self, mock_popen, mock_mkdir):
        """Start Anglerfish."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        af_sample_sheet = 'anglerfish_sample_sheet.csv'
        run = MinION(run_dir, None, af_sample_sheet)
        run.start_anglerfish()
        expected_parameters = ('anglerfish.py'
                          + ' --samplesheet anglerfish_sample_sheet.csv'
                          + ' --out_fastq data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2/anglerfish_output'
                          + ' --threads 2'
                          + ' --skip_demux'
                          + ' --skip_fastqc; echo $? > .exitcode_for_anglerfish')
        mock_popen.assert_called_once_with(expected_parameters, stdout=subprocess.PIPE, shell=True, cwd=run_dir)

    @mock.patch('taca.nanopore.minion.MinION._find_anglerfish_results')
    @mock.patch('taca.nanopore.minion.shutil.copyfile')
    def test_copy_results_for_lims(self, mock_copy, mock_results):
        """Copy Anglerfish results to lims."""
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        run = MinION(run_dir, None, None)
        anglerfish_results_path = 'anglerfish_output'
        anglerfish_results_file = os.path.join(run_dir, anglerfish_results_path, 'anglerfish_2020_09_23_141922', 'anglerfish_stats.txt')
        lims_results_file = 'some/dir/2020/anglerfish_stats_AAU644.txt'
        mock_results.return_value = anglerfish_results_file
        run.copy_results_for_lims()
        mock_copy.assert_called_once_with(anglerfish_results_file, lims_results_file)

    def test_find_anglerfish_results(self):
        """Locate Anglerfish results file."""
        anglerfish_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2/anglerfish_output'
        run_dir = 'data/nanopore_data/run4/done_demuxing/20200104_1412_MN19414_AAU644_68125dc2'
        run = MinION(run_dir, None, None)
        found_file = run._find_anglerfish_results()
        expected_file = os.path.join(anglerfish_dir, 'anglerfish_2020_09_23_141922', 'anglerfish_stats.txt')
        self.assertEqual(expected_file, found_file)
