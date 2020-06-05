#!/usr/bin/env python

import unittest
import mock

from taca.backup import backup
from taca.utils import config as conf

CONFIG = conf.load_yaml_config('data/taca_test_cfg_backup.yaml')


class TestRunVars(unittest.TestCase):
    """ Tests for TACA Backup variables class """

    def test_backup_variables(self):
        """ Set up backup variables """
        run_variables = backup.run_vars('data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX')
        self.assertEqual(run_variables.name, '190201_A00621_0032_BHHFCFDSXX')
        self.assertEqual(run_variables.zip, '190201_A00621_0032_BHHFCFDSXX.tar.gz')
        self.assertEqual(run_variables.key, '190201_A00621_0032_BHHFCFDSXX.key')
        self.assertEqual(run_variables.key_encrypted, '190201_A00621_0032_BHHFCFDSXX.key.gpg')
        self.assertEqual(run_variables.zip_encrypted, '190201_A00621_0032_BHHFCFDSXX.tar.gz.gpg')


class TestBackupUtils(unittest.TestCase):
    """ Tests for TACA Backup utils class """

    def test_fetch_config_info(self):
        """ Get backup info from config """
        config_info = backup.backup_utils('data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX')
        self.assertEqual(config_info.data_dirs, {'miseq': 'data/nas/miseq.lab'})
        self.assertEqual(config_info.archive_dirs, {'hiseq': 'blah', 'miseq': 'data/nas/miseq.lab/nosync'})
        self.assertEqual(config_info.keys_path, 'data/nas/run_keys')
        self.assertEqual(config_info.gpg_receiver, 'some.user')
        self.assertEqual(config_info.mail_recipients, 'some_user@some_email.com')
        self.assertEqual(config_info.check_demux, True)
        self.assertEqual(config_info.couch_info, {'url': 'url', 'username': 'username', 'password': 'pwd', 'port': 1234, 'xten_db': 'x_flowcells'})

    def test_collect_runs(self):
        """ Get backup runs from archive directories """
        backup_object = backup.backup_utils()
        backup_object.collect_runs(ext=".tar.gz", filter_by_ext=True)
        #import pdb; pdb.set_trace()
        run = backup_object.runs[0].name
        self.assertEqual(run, '200201_A00621_0032_BHHFCFDSXX')

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
        """ Check backup disk space """
        backup_object = backup.backup_utils()
        mock_sp.return_value = ["Filesystem   512-blocks      Used Available Capacity iused      ifree %iused  Mounted on\n/dev/disk1s1  976490576 100 813074776    15% 1086272 4881366608    0%   /System/Volumes/Data", None]
        path = 'data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX'
        run = '190201_A00621_0032_BHHFCFDSXX'
        with self.assertRaises(SystemExit):
            backup_object.avail_disk_space(path, run)

    @mock.patch('taca.backup.backup.sp.check_call')
    def test_file_in_pdc(self, mock_call):
        """ Check if files exist in PDC """
        mock_call.return_value = "Whatever"
        backup_object = backup.backup_utils()
        src_file = 'data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX/RTAComplete.txt'
        self.assertTrue(backup_object.file_in_pdc(src_file, silent=True))

    def test_get_run_type(self):
        """ Get run type from flowcell name """
        backup_object = backup.backup_utils()
        run_type = backup_object._get_run_type('190201_A00621_0032_BHHFCFDSXX')
        self.assertEqual(run_type, 'novaseq')

    def test_call_commands(self):
        """ Call expernal backup commands """
        backup_object = backup.backup_utils()
        got_output = backup_object._call_commands(cmd1="ls data/nas/miseq.lab", mail_failed=False, return_out=True)
        expected_output = (True, '190201_A00621_0032_BHHFCFDSXX\nnosync\n')
        self.assertEqual(got_output, expected_output)

    def test_check_status(self):
        """ Check subprocess status """
        backup_object = backup.backup_utils()
        cmd = 'ls'
        status_pass = 0
        err_msg = "Error"
        got_status_pass = backup_object._check_status(cmd, status_pass, err_msg, mail_failed=False)
        self.assertTrue(got_status_pass)
        status_fail = 1
        got_status_fail = backup_object._check_status(cmd, status_fail, err_msg, mail_failed=False)
        self.assertFalse(got_status_fail)

    @mock.patch('taca.backup.backup.os.remove')
    def test_clean_tmp_files(self, mock_remove):
        """Remove file if it exist"""
        backup_object = backup.backup_utils()
        files = ['data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX/RTAComplete.txt', 'data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX/missing_file.txt']
        backup_object._clean_tmp_files(files)
        mock_remove.assert_called_once_with('data/nas/miseq.lab/190201_A00621_0032_BHHFCFDSXX/RTAComplete.txt')
