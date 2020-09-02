#!/usr/bin/env python
import unittest
import mock
import crontab

from taca.server_status import server_status, cronjobs
from taca.utils import config

CONFIG = config.load_yaml_config('data/taca_test_cfg.yaml')
INITAL_TAB = """
# First Comment
0,30 * * * * firstcommand
"""

class TestServerStatus(unittest.TestCase):
    def test_get_nases_disk_space(self):
        """Get disk space for disk specified in config file."""
        got_disk_space = server_status.get_nases_disk_space()
        assert got_disk_space

    def test_parse_output_valid_case(self):
        """Parse valid disk space output."""
        valid_disk_space = 'Filesystem     Size   Used  Avail Capacity iused      ifree %iused  Mounted on \
        /dev/disk1s1  466Gi   59Gi  393Gi    14% 1062712 4881390168    0%   /System/Volumes/Data'
        expected_result = {'disk_size': '14%',
                           'mounted_on': '/System/Volumes/Data',
                           'available_percentage': '100%',
                           'space_used': '1062712',
                           'used_percentage': '0%',
                           'filesystem': '393Gi',
                           'space_available': '4881390168'}
        got_result = server_status._parse_output(valid_disk_space)
        self.assertEqual(expected_result, got_result)

    def test_parse_output_invalid_case(self):
        """Parse invalid disk space output."""
        invalid_disk_space = ''
        expected_invalid_result = {
            'disk_size': 'NaN',
            'space_used': 'NaN',
            'space_available': 'NaN',
            'used_percentage': 'NaN',
            'available_percentage': 'NaN',
            'mounted_on': 'NaN',
            'filesystem': 'NaN'
        }
        invalid_result = server_status._parse_output(invalid_disk_space)
        self.assertEqual(expected_invalid_result, invalid_result)

    @mock.patch('taca.server_status.server_status.statusdb')
    def test_update_status_db(self, mock_couchdb):
        """Update statusdb."""
        disk_space = {'localhost': {'disk_size': '14%', 'mounted_on': '/System/Volumes/Data', 'available_percentage': '100%', 'space_used': '1061701', 'used_percentage': '0%', 'filesystem': '393Gi', 'space_available': '4881391179'}}
        server_status.update_status_db(disk_space, server_type='nas')


class TestCronjobs(unittest.TestCase):
    @mock.patch('taca.server_status.cronjobs.CronTab')
    @mock.patch('taca.server_status.cronjobs.getpass.getuser')
    def test_parse_crontab(self, mock_getpass, mock_crontab):
        """Parse crontab."""
        mock_crontab.return_value = crontab.CronTab(tab=INITAL_TAB)
        mock_getpass.return_value = 'test_user'
        expected_crontab = {'test_user':
                            [{'Comment': u'First Comment',
                              'Day of month': '*',
                              'Command': u'firstcommand',
                              'Hour': '*',
                              'Day of week': '*',
                              'Enabled': True,
                              'Special syntax': '',
                              'Minute': '0,30',
                              'Month': '*'}]
        }

        got_crontab = cronjobs._parse_crontab()
        self.assertEqual(expected_crontab, got_crontab)

    @mock.patch('taca.server_status.cronjobs.statusdb')
    @mock.patch('taca.server_status.cronjobs.logging')
    @mock.patch('taca.server_status.cronjobs.platform')
    @mock.patch('taca.server_status.cronjobs._parse_crontab')
    def test_update_cronjob_db(self, mock_parser, mock_platform, mock_logging, mock_statusdb):
        """Update couchdb with cronjobs."""
        mock_parser.return_value = {'test_user':
                            [{'Comment': u'First Comment',
                              'Day of month': '*',
                              'Command': u'firstcommand',
                              'Hour': '*',
                              'Day of week': '*',
                              'Enabled': True,
                              'Special syntax': '',
                              'Minute': '0,30',
                              'Month': '*'}]
        }
        mock_platform.node.return_value = 'server.name'
        cronjobs.update_cronjob_db()
        calls = [mock.call.info('Connecting to database: url'),
                 mock.call.warning('Document has not been created/updated')]
        mock_logging.assert_has_calls(calls)
