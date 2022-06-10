#!/usr/bin/env python
import unittest
import mock
import tempfile
import os
import shutil

from taca.backup import backup
from taca.utils import config as conf

CONFIG = conf.load_yaml_config('data/taca_test_cfg_backup.yaml')


class TestRunVars(unittest.TestCase):
    """Tests for TACA Backup variables class."""

    def test_backup_variables(self):
        """Set up backup variables."""
        run_variables = backup.run_vars('data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX')
        self.assertEqual(run_variables.name, '190201_A00621_0032_BHHFCFDSXX')
        self.assertEqual(run_variables.zip, '190201_A00621_0032_BHHFCFDSXX.tar.gz')
        self.assertEqual(run_variables.key, '190201_A00621_0032_BHHFCFDSXX.key')
        self.assertEqual(run_variables.key_encrypted, '190201_A00621_0032_BHHFCFDSXX.key.gpg')
        self.assertEqual(run_variables.zip_encrypted, '190201_A00621_0032_BHHFCFDSXX.tar.gz.gpg')


class TestBackupUtils(unittest.TestCase):
    """Tests for TACA Backup utils class."""

    def test_fetch_config_info(self):
        """Get backup info from config."""
        config_info = backup.backup_utils('data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX')
        self.assertEqual(config_info.data_dirs, {'miseq': 'data/nas/miseq.lab'})
        self.assertEqual(config_info.archive_dirs, {'hiseq': 'blah', 'miseq': 'data/nas/miseq.lab/nosync'})
        self.assertEqual(config_info.keys_path, 'data/nas/run_keys')
        self.assertEqual(config_info.gpg_receiver, 'some.user')
        self.assertEqual(config_info.mail_recipients, 'some_user@some_email.com')
        self.assertEqual(config_info.check_demux, True)
        self.assertEqual(config_info.couch_info, {'url': 'url', 'username': 'username', 'password': 'pwd', 'xten_db': 'x_flowcells'})

    def test_collect_runs(self):
        """Get backup runs from archive directories."""
        backup_object = backup.backup_utils()
        backup_object.collect_runs(ext='.tar.gz', filter_by_ext=True)
        run = backup_object.runs[0].name
        self.assertEqual(run, '200201_A00621_0032_BHHFCFDSXY')

    def test_collect_runs_specific_run(self):
        """Collect only specific run."""
        backup_object = backup.backup_utils(run='data/nas/miseq.lab/nosync/200201_A00621_0032_BHHFCFDSXX')
        backup_object.collect_runs()
        run = backup_object.runs[0].name
        self.assertEqual(run, '200201_A00621_0032_BHHFCFDSXX')

        missing_object = backup.backup_utils(run='some/missing/path/run')
        with self.assertRaises(SystemExit):
            missing_object.collect_runs()

    @mock.patch('taca.backup.backup.sp.Popen.communicate')
    @mock.patch('taca.backup.backup.misc')
    def test_avail_disk_space(self, mock_misc, mock_sp):
        """Check backup disk space."""
        backup_object = backup.backup_utils()
        mock_sp.return_value = ['Filesystem   512-blocks      Used Available Capacity iused      ifree %iused  Mounted on\n/dev/disk1s1  976490576 100 813074776    15% 1086272 4881366608    0%   /System/Volumes/Data', None]
        path = 'data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX'
        run = '190201_A00621_0032_BHHFCFDSXX'
        with self.assertRaises(SystemExit):
            backup_object.avail_disk_space(path, run)

    @mock.patch('taca.backup.backup.sp.check_call')
    def test_file_in_pdc(self, mock_call):
        """Check if files exist in PDC."""
        mock_call.return_value = 'Whatever'
        backup_object = backup.backup_utils()
        src_file = 'data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX/RTAComplete.txt'
        self.assertTrue(backup_object.file_in_pdc(src_file, silent=True))

    def test_get_run_type(self):
        """Get run types from flowcell names."""
        backup_object = backup.backup_utils()
        novaseq_run = backup_object._get_run_type('190201_A00621_0032_BHHFCFDSXX')
        self.assertEqual(novaseq_run, 'novaseq')
        hiseqx_run = backup_object._get_run_type('190711_ST-E00266_0356_AH2L32CCX2')
        self.assertEqual(hiseqx_run, 'hiseqx')
        miseq_run = backup_object._get_run_type('200604_M01320_0329_000000000-J668J')
        self.assertEqual(miseq_run, 'miseq')
        hiseq_run = backup_object._get_run_type('190628_D00415_0465_BH2HVYBCX3')
        self.assertEqual(hiseq_run, 'hiseq')
        nextseq_run = backup_object._get_run_type('200602_NS500688_0656_AHGCKWBGXF')
        self.assertEqual(nextseq_run, 'nextseq')

    def test_call_commands(self):
        """Call expernal backup command."""
        backup_object = backup.backup_utils()
        got_output = backup_object._call_commands(cmd1='ls data/nas/miseq.lab', mail_failed=False, return_out=True)
        expected_output = (True, b'190201_A00621_0032_BHHFCFDSXX\nnosync\n')
        self.assertEqual(got_output, expected_output)

    def test_call_commands_double(self):
        """Call external backup command, given two commands."""
        backup_object = backup.backup_utils()
        tmp_dir = os.path.join(tempfile.mkdtemp(), 'tmp')
        tmp_file = os.path.join(tmp_dir, 'output.out')
        os.makedirs(tmp_dir)
        cmd1 = 'ls data/nas/miseq.lab'
        cmd2 = 'ls data/nas/miseq.lab'
        got_output = backup_object._call_commands(cmd1, cmd2, out_file=tmp_file, mail_failed=False)
        self.assertTrue(os.path.isfile(tmp_file))
        shutil.rmtree(tmp_dir)

    def test_check_status(self):
        """Check subprocess status."""
        backup_object = backup.backup_utils()
        cmd = 'ls'
        status_pass = 0
        err_msg = 'Error'
        got_status_pass = backup_object._check_status(cmd, status_pass, err_msg, mail_failed=False)
        self.assertTrue(got_status_pass)
        status_fail = 1
        got_status_fail = backup_object._check_status(cmd, status_fail, err_msg, mail_failed=False)
        self.assertFalse(got_status_fail)

    @mock.patch('taca.backup.backup.os.remove')
    def test_clean_tmp_files(self, mock_remove):
        """Remove file if it exist."""
        backup_object = backup.backup_utils()
        files = ['data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX/RTAComplete.txt', 'data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX/missing_file.txt']
        backup_object._clean_tmp_files(files)
        mock_remove.assert_called_once_with('data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX/RTAComplete.txt')

    @mock.patch('taca.backup.backup.statusdb', autospec=True)
    @mock.patch('taca.backup.backup.logger')
    def test_log_pdc_statusdb(self, mock_logger, mock_couch):
        """Update statusdb if transfer was successful."""
        backup_object = backup.backup_utils()
        run = '190201_A00621_0032_BHHFCFDSXX'
        backup_object._log_pdc_statusdb(run)
        mock_logger.warn.assert_called_once()

    @mock.patch('taca.backup.backup.backup_utils._call_commands', return_value=True)
    @mock.patch('taca.backup.backup.shutil')
    @mock.patch('taca.backup.backup.backup_utils._clean_tmp_files')
    @mock.patch('taca.backup.backup.backup_utils.avail_disk_space')
    def test_encrypt_runs(self, mock_space, mock_clean, mock_shutil, mock_command):
        """Encrypt found runs."""
        backup_object = backup.backup_utils(run='data/nas/miseq.lab/nosync/200201_A00621_0032_BHHFCFDSXX')
        run = 'data/nas/miseq.lab/nosync/190201_A00621_0032_BHHFCFDSXX'
        force = True
        backup_object.encrypt_runs(run, force)
        mock_clean.assert_called_once()
        os.remove('data/nas/miseq.lab/nosync/190201_A00621_0032_BHHFCFDSXX.encrypting')

    @mock.patch('taca.backup.backup.logger.error')
    def test_pdc_put(self, mock_logger):
        """Put runs on PDC."""
        backup_object = backup.backup_utils(run='data/nas/miseq.lab/nosync/200201_A00621_0032_BHHFCFDSXX')
        run = 'data/nas/miseq.lab/nosync/190201_A00621_0032_BHHFCFDSXX'
        backup_object.pdc_put(run)
        mock_logger.assert_called_once()
