#!/usr/bin/env python

import os
import shutil
import tempfile
import unittest
import mock
from datetime import datetime

from taca.cleanup import cleanup
from taca.utils import config as conf

CONFIG = conf.load_yaml_config('data/taca_test_cfg_cleanup.yaml')


class TestCleanup(unittest.TestCase):
    """Tests for TACA Cleanup module."""

    @mock.patch('taca.cleanup.cleanup.statusdb')
    @mock.patch('taca.cleanup.cleanup.get_closed_proj_info')
    @mock.patch('taca.cleanup.cleanup.misc.query_yes_no')
    @mock.patch('taca.cleanup.cleanup._remove_files')
    @mock.patch('taca.cleanup.cleanup._touch_cleaned')
    def test_cleanup_miarka(self, mock_touch, mock_rm, mock_query, mock_info,  mock_statusdb):
        """Locate and move old data on Miarka."""
        mock_info.return_value = {'closed_date': '2019-04-07',
                                  'bioinfo_responsible': 'O.B. One',
                                  'pid': 'P1234',
                                  'name': 'N.Owens_19_01',
                                  'closed_days': 5}
        mock_query.return_value = True
        mock_rm.return_value = True
        days_fastq = 1
        days_analysis = 1
        only_fastq = False
        only_analysis = False
        clean_undetermined = False
        status_db_config = 'data/taca_test_cfg_cleanup.yaml'
        exclude_projects = False
        list_only = False
        date = '2016-01-31'
        calls = [mock.call('data/miarka/incoming/190201_A00621_0032_BHHFCFDSXX/Demultiplexing/N.Owens_19_01'),
                 mock.call('../../nobackup/NGI/ANALYSIS/P1234')]
        cleanup.cleanup_miarka(days_fastq, days_analysis, only_fastq, only_analysis, clean_undetermined, status_db_config, exclude_projects, list_only, date, dry_run=False)
        mock_touch.assert_has_calls(calls)

    def test_get_closed_proj_info(self):
        """Return a dict if project is closed."""
        pid = 'P1234'
        pdoc = {'close_date': '2019-04-07',
                'project_name': 'A.Name_19_01',
                'project_id': 'P1234',
                'project_summary': {'bioinfo_responsible': 'O.B. One'}}
        tdate = datetime.strptime('2019-04-08', '%Y-%m-%d')
        got_data = cleanup.get_closed_proj_info(pid, pdoc, tdate)
        expected_data = {'closed_date': '2019-04-07',
                         'bioinfo_responsible': 'O.B. One',
                         'pid': 'P1234',
                         'name': 'A.Name_19_01',
                         'closed_days': 1}
        self.assertEqual(got_data, expected_data)

    def test_collect_analysis_data_miarka(self):
        """Get analysis data on Miarka."""
        pid = 'P1234'
        analysis_root = 'data/test_data/analysis'
        file_list, size = cleanup.collect_analysis_data_miarka(pid, analysis_root, files_ext_to_remove={})
        self.assertEqual(file_list, 'cleaned')

    def test_collect_fastq_data_miarka(self):
        """Collect removed files."""
        fc_root = 'data/test_data/190201_A00621_0032_BHHFCFDSXX'
        fc_proj_src = 'N.Owens_19_01'
        file_list, size = cleanup.collect_fastq_data_miarka(fc_root, fc_proj_src)
        expected_data = {'flowcells':
                         {'190201_A00621_0032_BHHFCFDSXX':
                          {'proj_root': 'data/test_data/190201_A00621_0032_BHHFCFDSXX/N.Owens_19_01',
                           'fq_files': ['data/test_data/190201_A00621_0032_BHHFCFDSXX/N.Owens_19_01/sample1.fastq.gz',
                                        'data/test_data/190201_A00621_0032_BHHFCFDSXX/N.Owens_19_01/sample2.fastq.gz']}}}
        self.assertEqual(file_list, expected_data)
        self.assertEqual(size, 0)

    def test_collect_files_by_ext(self):
        """Return found paths."""
        path = 'data/test_data'
        ext = ['*.txt']
        found_files = cleanup.collect_files_by_ext(path, ext)
        expected_files = ['data/test_data/nosync/190201_A00621_0032_BHHFCFDSXY/RTAComplete.txt',
                          'data/test_data/190201_A00621_0032_BHHFCFDSXX/RTAComplete.txt']
        self.assertEqual(sorted(found_files), sorted(expected_files))

    def test_get_proj_meta_info(self):
        """Get project metadata."""
        info = {'name': 'Nobody Owens',
                'pid': 'P1234',
                'bioinfo_responsible': 'O.B. One',
                'closed_days': 1,
                'closed_date': '2020-04-07',
                'fastq_size': 1001}
        days_fastq = ''
        got_data = cleanup.get_proj_meta_info(info, days_fastq)
        expected_data = '''
Project overview: Nobody Owens
Project ID: P1234
Bioinfo Responsible: O.B. One
Closed for (days): 1
Closed from (date): 2020-04-07
Project analysis: No analysis directory
Estimated data size: ~2kb
'''
        self.assertEqual(got_data, expected_data)

    def test_get_files_size_text(self):
        """Format file size string."""
        plist = {'P1': {'fastq_size': 1001, 'analysis_size': 1000000},
                 'P2': {'fastq_size': 1001, 'analysis_size': 1000000}}
        got_data = cleanup.get_files_size_text(plist)
        expected_data = '(~~2kb fastq data and ~~2mb analysis data) '
        self.assertEqual(got_data, expected_data)

    def test_def_get_size_unit(self):
        """Convert size."""
        #function broken if size < 1000
        size = 1001
        self.assertEqual(cleanup._def_get_size_unit(size), '~1kb')
        size *= 1000
        self.assertEqual(cleanup._def_get_size_unit(size), '~1mb')
        size *= 1000
        self.assertEqual(cleanup._def_get_size_unit(size), '~1gb')
        size *= 1000
        self.assertEqual(cleanup._def_get_size_unit(size), '~1tb')

    @mock.patch('taca.cleanup.cleanup.os.remove')
    def test_remove_files(self, mock_remove):
        """Remove files in given list."""
        files = ['file1', 'file2']
        cleanup._remove_files(files)
        calls = [mock.call('file1'), mock.call('file2')]
        mock_remove.assert_has_calls(calls)

    def test_touch_cleaned(self):
        """Create empty file in specified dir."""
        tmp_dir = os.path.join(tempfile.mkdtemp(), 'tmp')
        os.makedirs(tmp_dir)
        cleanup._touch_cleaned(tmp_dir)
        expected_file = os.path.join(tmp_dir, 'cleaned')
        self.assertTrue(os.path.exists(expected_file))
        shutil.rmtree(tmp_dir)
