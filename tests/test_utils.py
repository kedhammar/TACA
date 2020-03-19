""" Unit tests for the utils helper functions """

import hashlib
import mock
import os
import shutil
import subprocess
import tempfile
import unittest
import time
import couchdb
from taca.utils import misc, filesystem, transfer, config, bioinfo_tab


class TestMisc(unittest.TestCase):
    """ Test class for the misc functions """

    @classmethod
    def setUpClass(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_taca_misc")
        self.hashfile = os.path.join(self.rootdir,'test_hashfile')
        with open(self.hashfile,'w') as fh:
            fh.write("This is some contents\n")
        self.hashfile_digests = {
            'SHA256':
                '4f075ae76b480bb0200dab01cd304f4045e04cd2b73e88b89549e5ac1627f222',
            'MD5':
                'c8498fc299bc3e22690045f1b62ce4e9',
            'SHA1':
                '098fb272dfdae2ea1ba57c795dd325fa70e3c3fb'}

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.rootdir)

    # Test generator for different hashing algorithms
    def test_hashfile(self):
        for alg,obj in self.hashfile_digests.items():
            yield self.check_hash, alg, obj

    def test_hashfile_dir(self):
        """Hash digest for a directory should be None"""
        assert misc.hashfile(self.rootdir) is None

    def test_multiple_hashfile_calls(self):
        """ Ensure that the hasher object is cleared between subsequent calls
        """
        assert misc.hashfile(self.hashfile,hasher='sha1') == misc.hashfile(self.hashfile,'sha1')

    def check_hash(self, alg, exp):
        assert misc.hashfile(self.hashfile,hasher=alg) == exp

    @mock.patch('taca.utils.misc.smtplib.SMTP')
    def test_send_mail(self, mock_smtplib):
        """ Test send email """
        mock_smtplib.sendmail()
        mock_smtplib.quit()
        assert misc.send_mail("subject", "content", "receiver") is None
        with self.assertRaises(SystemExit):
            misc.send_mail("subject", "content", None)

    def test_call_external_command_pass(self):
        """ Call external command """
        new_file = os.path.join(self.rootdir, "test_call_external")
        command = "touch " + new_file
        misc.call_external_command(command)
        assert os.path.isfile(new_file)

    def test_call_external_command_fail(self):
        """ Call external command should handle error """
        command = "ls -Q"
        with self.assertRaises(subprocess.CalledProcessError):
            misc.call_external_command(command)

    def test_call_external_command_detached(self):
        """ Call external command detached"""
        new_file = os.path.join(self.rootdir, "test_call_external")
        command = "touch " + new_file
        misc.call_external_command_detached(command)
        time.sleep(0.1)
        assert os.path.isfile(new_file)

    def test_to_seconds(self):
        """ Transform days and hours to seconds """
        with self.assertRaises(SystemExit):
            misc.to_seconds(days=1, hours=1)
        with self.assertRaises(SystemExit):
            misc.to_seconds()
        self.assertEqual(misc.to_seconds(days=1), 86400)
        self.assertEqual(misc.to_seconds(hours=1), 3600)

    @mock.patch('taca.utils.misc.raw_input', return_value='yes')
    def test_query_yes_no_true(self, mock_raw_input):
        """Return True from answer yes """
        response = misc.query_yes_no("Some question")
        self.assertTrue(response)

    @mock.patch('taca.utils.misc.raw_input', return_value='no')
    def test_query_yes_no_false(self, mock_raw_input):
        """Return False from answer no """
        response = misc.query_yes_no("Some question")
        self.assertFalse(response)

    def test_return_unique(self):
        """ Return unique items in a list """
        input_list = ['a', 'b', 'a', 'c']
        returned_list = misc.return_unique(input_list)
        expected_list = ['a', 'b', 'c']
        self.assertEqual(returned_list, expected_list)


class TestFilesystem(unittest.TestCase):
    """ Test class for the filesystem functions """

    def setUp(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_taca_filesystem")

    def tearDown(self):
        shutil.rmtree(self.rootdir)

    def test_crete_folder_non_existing(self):
        """ Ensure that a non-existing folder is created """
        target_folder = os.path.join(self.rootdir,"target-non-existing")
        self.assertTrue(
            filesystem.create_folder(target_folder),
            "A non-existing target folder could not be created")
        self.assertTrue(
            os.path.exists(target_folder),
            "A non-existing target folder was not created \
            but method returned True"
        )

    def test_crete_folder_existing(self):
        """ Ensure that an existing folder is detected """
        self.assertTrue(
            filesystem.create_folder(self.rootdir),
            "A pre-existing target folder was not detected")

    def test_crete_folder_parent_non_existing(self):
        """ Ensure that a non-existing parent folder is created """
        target_folder = os.path.join(
            self.rootdir,
            "parent-non-existing",
            "target-non-existing")
        self.assertTrue(
            filesystem.create_folder(target_folder),
            "A non-existing parent and target folder could not be created")
        self.assertTrue(
            os.path.exists(target_folder),
            "A non-existing parent folder was not created \
            but method returned True"
        )

    def test_crete_folder_exception(self):
        """ Ensure that create_folder handles thrown exceptions gracefully """
        with mock.patch.object(filesystem.os, 'makedirs', side_effect=OSError):
            self.assertFalse(
                filesystem.create_folder(
                    os.path.join(self.rootdir,"target-non-existing")),
                "A raised exception was not handled properly")

    def test_chdir(self):
        """ Ensure start dir and end dir are the same """
        initial_dir = os.getcwd()
        filesystem.chdir(self.rootdir)
        final_dir = os.getcwd()
        self.assertEqual(initial_dir, final_dir)


class TestTransferAgent(unittest.TestCase):
    """ Test class for the TransferAgent class """

    @classmethod
    def setUpClass(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_taca_transfer_src")
        self.testfile = tempfile.mkstemp(dir=self.rootdir)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.rootdir)

    def setUp(self):
        self.destdir = tempfile.mkdtemp(prefix="test_taca_transfer_dest")
        self.agent = transfer.TransferAgent(
            src_path=self.rootdir,
            dest_path=self.destdir)

    def tearDown(self):
        shutil.rmtree(self.destdir)

    def test_transfer_validate_src_path(self):
        """ src_path should validate properly """
        self.agent.validate_src_path()
        self.agent.src_path = None
        with self.assertRaises(transfer.TransferError):
            self.agent.validate_src_path()
        self.agent.src_path = os.path.join(
            self.rootdir,
            "this-file-does-not-exist")
        with self.assertRaises(transfer.TransferError):
            self.agent.validate_src_path()

    def test_transfer_validate_dest_path(self):
        """ dest_path should validate properly """
        self.agent.validate_dest_path()
        self.agent.dest_path = None
        with self.assertRaises(transfer.TransferError):
            self.agent.validate_dest_path()

    def test_transfer_transfer(self):
        """ do_transfer in superclass should raise exception if called """
        with self.assertRaises(NotImplementedError):
            self.agent.transfer()

    def test_transfer_validate_transfer(self):
        """ validate_transfer in superclass should raise exception if called """
        with self.assertRaises(NotImplementedError):
            self.agent.validate_transfer()


class TestSymlinkAgent(unittest.TestCase):
    """ Test class for the SymlinkAgent class """

    @classmethod
    def setUpClass(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_taca_symlink_src")
        path = self.rootdir
        for n in xrange(3):
            open(os.path.join(path,"file{}".format(n)),'w').close()
            path = os.path.join(path,"folder{}".format(n))
            os.mkdir(path)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.rootdir)

    def setUp(self):
        self.targetdir = tempfile.mkdtemp(
            prefix="test_taca_filesystem_symlink_dest")

    def tearDown(self):
        shutil.rmtree(self.targetdir)

    def test_symlink_validate_transfer(self):
        """ verify that the dest_path was created """
        src = os.path.join(self.rootdir,"file0")
        dst = os.path.join(self.targetdir,"file0")
        os.symlink(src,dst)
        self.assertTrue(transfer.SymlinkAgent(src,dst).validate_transfer())

    def test_symlink_file_top_folder(self):
        """ Symlink a single file in the top folder """
        src = os.path.join(self.rootdir,"file0")
        target = os.path.join(self.targetdir,os.path.basename(src))
        self.assertTrue(transfer.SymlinkAgent(src,target).transfer())

    def test_symlink_file_make_dir(self):
        """ Symlink a single file into a non-existing folder """
        src = os.path.join(self.rootdir,"folder0","folder1","file2")
        target = os.path.join(
            self.targetdir,
            "these","folders","should","be","created")
        self.assertTrue(transfer.SymlinkAgent(src,target).transfer())

    def test_symlink_file_overwrite(self):
        """ Replace an existing file with overwrite """
        src = os.path.join(self.rootdir,"file0")
        target = os.path.join(self.targetdir,os.path.basename(src))
        open(target,'w').close()
        self.assertTrue(transfer.SymlinkAgent(src,target).transfer())

    def test_symlink_file_not_overwrite(self):
        """ Don't replace an existing file without overwrite """
        src = os.path.join(self.rootdir,"file0")
        target = os.path.join(self.targetdir,os.path.basename(src))
        open(target,'w').close()
        self.assertFalse(
            transfer.SymlinkAgent(src,target,overwrite=False).transfer())

    def test_symlink_file_broken(self):
        """ Don't create a broken symlink """
        src = os.path.join(self.rootdir,"non-existing-file")
        target = os.path.join(self.targetdir,os.path.basename(src))
        with self.assertRaises(transfer.TransferError):
            transfer.SymlinkAgent(src,target).transfer()

    def test_symlink_file_unlink_fail(self):
        """ Failing to remove existing file should raise SymlinkError """
        src = self.rootdir
        target = os.path.join(self.targetdir,"target-file")
        open(target,'w').close()
        with mock.patch.object(
            transfer.os,
            'unlink',
            side_effect=OSError("Mocked error")):
            with self.assertRaises(transfer.SymlinkError):
                transfer.SymlinkAgent(src,target).transfer()

    def test_symlink_folder_top_folder(self):
        """ Symlinking a top-level folder """
        src = os.path.join(self.rootdir,"folder0")
        target = os.path.join(self.targetdir,os.path.basename(src))
        self.assertTrue(transfer.SymlinkAgent(src,target).transfer())

    def test_symlink_folder_overwrite(self):
        """ Replace an existing folder with overwrite """
        src = os.path.join(self.rootdir,"folder0")
        target = os.path.join(self.targetdir,os.path.basename(src))
        shutil.copytree(src,target)
        self.assertTrue(transfer.SymlinkAgent(src,target).transfer())

    def test_symlink_folder_mount_point(self):
        """ Don't overwrite a mount point """
        src = os.path.join(self.rootdir)
        target = os.path.join(self.targetdir)
        with mock.patch.object(transfer.os.path,'ismount',return_value=True):
            with self.assertRaises(transfer.SymlinkError):
                transfer.SymlinkAgent(src,target).transfer()

    def test_symlink_folder_not_overwrite(self):
        """ Don't overwrite an existing path that is neither a mount point,
            file, link or directory
        """
        src = os.path.join(self.rootdir)
        target = os.path.join(self.targetdir)
        with mock.patch('taca.utils.transfer.os.path') as mockobj:
            mockobj.ismount.return_value = False
            mockobj.isfile.return_value = False
            mockobj.islink.return_value = False
            mockobj.isdir.return_value = False
            with self.assertRaises(transfer.SymlinkError):
                transfer.SymlinkAgent(src,target).transfer()

    def test_symlink_folder_parent_error(self):
        """ Failing to create parent folder structure should raise SymlinkError
        """
        src = self.rootdir
        target = os.path.join(self.targetdir,"non-existing-folder","target-file")
        with mock.patch.object(transfer,'create_folder',return_value=False):
            with self.assertRaises(transfer.SymlinkError):
                transfer.SymlinkAgent(src,target).transfer()

    def test_symlink_folder_rmtree_fail(self):
        """ Failing to remove existing folder should raise SymlinkError
        """
        src = self.rootdir
        target = self.targetdir
        with mock.patch.object(
            transfer.shutil,
            'rmtree',
            side_effect=OSError("Mocked error")):
            with self.assertRaises(transfer.SymlinkError):
                transfer.SymlinkAgent(src,target).transfer()

    def test_symlink_folder_symlink_error(self):
        """ Failing to create symlink should raise SymlinkError
        """
        src = self.rootdir
        target = os.path.join(self.targetdir,os.path.basename(src))
        with mock.patch.object(
            transfer.os,
            'symlink',
            side_effect=OSError("Mocked error")):
            with self.assertRaises(transfer.SymlinkError):
                transfer.SymlinkAgent(src,target).transfer()

    def test_symlink_folder_unexpected(self):
        """ An unexpected exception should propagate upwards
        """
        src = self.rootdir
        target = self.targetdir
        with mock.patch.object(
            transfer.os.path,
            'exists',
            side_effect=Exception("Mocked error")):
            with self.assertRaises(Exception):
                transfer.SymlinkAgent(src,target).transfer()


class TestRsyncAgent(unittest.TestCase):
    """ Test class for the RsyncAgent class """

    @classmethod
    def setUpClass(cls):
        cls.rootdir = tempfile.mkdtemp(prefix="test_taca_transfer_src")
        (fh, cls.testfile) = tempfile.mkstemp(
            prefix="test_taca_transfer_file")
        os.write(fh,"this is some content")
        os.close(fh)
        open(os.path.join(cls.rootdir,"file0"),'w').close()
        f = os.path.join(cls.rootdir,"folder0")
        os.mkdir(f)
        open(os.path.join(f,"file1"),'w').close()

        # create a digest file
        def _write_digest(rootdir, fhandle, fpath):
            fhandle.write("{}  {}\n".format(misc.hashfile(fpath), os.path.relpath(fpath, rootdir)))

        cls.digestfile = os.path.join(cls.rootdir, 'digestfile.sha1')
        with open(cls.digestfile, 'w') as digesth:
            map(lambda x:
                map(lambda y: _write_digest(cls.rootdir, digesth, os.path.join(x[0], y)),
                    filter(lambda z: os.path.join(x[0], z) != cls.digestfile, x[2])),
                os.walk(cls.rootdir))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.rootdir)
        os.unlink(cls.testfile)

    def setUp(self):
        self.destdir = tempfile.mkdtemp(prefix="test_taca_transfer_dest")
        self.agent = transfer.RsyncAgent(
            self.rootdir,
            dest_path=self.destdir,
            validate=False)

    def tearDown(self):
        shutil.rmtree(self.destdir)

    def test_init(self):
        """ test initiation of agent instance """
        # initiate with some values
        args = ["arg1"]
        kwargs = {"dest_path": "arg2",
                  "remote_host": "arg3",
                  "remote_user": "arg4",
                  "validate": True,
                  "digestfile": "arg5"}
        agent = transfer.RsyncAgent(*args, **kwargs)
        # assert that the initialized values are what is set on the instance
        self.assertEqual(getattr(agent,"src_path"),args[0])
        for attribute, value in kwargs.items():
            self.assertEqual(getattr(agent, attribute), value)
        self.assertEqual(agent.cmdopts, agent.DEFAULT_OPTS)

    def test_rsync_validate_transfer(self):
        """ validate_transfer """
        # validation on remote hosts are not supported
        self.agent.remote_host = "not None"
        with self.assertRaises(NotImplementedError):
            self.agent.validate_transfer()
        # validation without a digestfile throws an exception
        self.agent.remote_host = None
        with self.assertRaises(transfer.RsyncValidationError):
            self.agent.validate_transfer()
        # validation with a valid digestfile should return true
        self.agent.digestfile = self.digestfile
        self.assertTrue(self.agent.validate_transfer(), "validation with a valid digestfile should return true")
        # modifying the contents of the digestfile should make validation fail
        with open(self.digestfile, 'a') as fh:
            fh.write("randomdigeststring  this-file-does-not-exist")
        self.assertFalse(self.agent.validate_transfer(), "validation with an invalid digestfile should return false")

    def test_rsync_validate_dest_path(self):
        """ Destination path should be properly checked """
        try:
            self.agent.validate_dest_path()
        except transfer.TransferError as e:
            self.fail("a proper path raised an exception: {}".format(e))
        self.agent.remote_host = None
        self.agent.dest_path = None
        with self.assertRaises(transfer.TransferError):
            self.agent.validate_dest_path()
        self.agent.remote_user = "user"
        self.agent.dest_path = self.destdir
        with self.assertRaises(transfer.TransferError):
            self.agent.validate_dest_path()

    def test_rsync_agent_dest_paths_constructed(self):
        """ Destination path should be properly constructed """
        self.assertEqual(
            self.destdir,
            self.agent.remote_path(),
            "Destination path was not correct for empty remote user " \
            "and empty destination host")
        self.agent.remote_host = "localhost"
        self.assertEqual(
            "localhost:{}".format(self.destdir),
            self.agent.remote_path(),
            "Destination path was not correct for empty remote user")
        self.agent.remote_user = "user"
        self.assertEqual(
            "user@localhost:{}".format(self.destdir),
            self.agent.remote_path(),
            "Destination path was not correct for non-empty remote user")
        self.agent.dest_path = None
        self.assertEqual(
            "user@localhost:",
            self.agent.remote_path(),
            "Destination path was not correct for empty destination path")

    def test_rsync_agent_propagate_error(self):
        """ An error thrown by the rsync subprocess should be wrapped and
            propagated """
        with mock.patch.object(
            transfer.subprocess,'check_call',
                side_effect=subprocess.CalledProcessError(
                    cmd="mocked subprocess",
                    returncode=-1)):
            with self.assertRaises(transfer.RsyncError):
                self.agent.transfer()

    def test_rsync_agent_file(self):
        """ rsync transfer of a single file """
        self.agent.src_path = os.path.join(self.rootdir,"file0")
        self.assertTrue(
            self.agent.transfer(),
            "transfer a single file failed")
        self.assertTrue(
            self.validate_files(
                self.agent.src_path,
                os.path.join(
                    self.destdir,
                    os.path.basename(self.agent.src_path))),
            "test file was not properly transferred")

    def test_rsync_agent_dir(self):
        """ rsync transfer of a folder """
        self.agent.src_path = os.path.join(self.rootdir,"folder0")
        self.assertTrue(
            self.agent.transfer(),
            "transfer a folder failed")
        self.assertTrue(
            self.validate_folders(
                self.agent.src_path,
                os.path.join(
                    self.destdir,
                    os.path.basename(self.agent.src_path))),
            "folder was not properly transferred")

    def test_rsync_agent_symlink(self):
        """ rsync should be able to resolve symlinks """
        self.agent.src_path = os.path.join(self.rootdir,"folder0")
        os.symlink(self.testfile,os.path.join(self.agent.src_path,"link1"))
        self.agent.cmdopts = {'-a': None, '--copy-links': None}
        self.assertTrue(
            self.agent.transfer(),
            "transfer a folder containing a symlink failed")
        self.assertEqual(
            misc.hashfile(self.testfile,hasher='sha1'),
            misc.hashfile(
                os.path.join(self.destdir,"folder0","link1"),
                hasher='sha1'),
            "symlink was not properly transferred")

    def validate_folders(self,src,dst):
        for root, dirs, files in os.walk(src):
            for file in files:
                s = os.path.join(root,file)
                d = os.path.join(dst,os.path.relpath(s,src))
                if not self.validate_files(s,d):
                    return False
        return True

    def validate_files(self,src,dst):
        return os.path.exists(src) and \
            os.path.isfile(src) and \
            os.path.exists(dst) and \
            os.path.isfile(dst) and \
            misc.hashfile(src) == misc.hashfile(dst)

class TestConfig(unittest.TestCase):

    def test_load_yaml_config_pass(self):
        """ Load a yaml config file """
        got_config_data = config.load_yaml_config('data/taca_test_cfg_minimal.yaml')
        expexted_config_data = {'statusdb':
                                {'url': 'url',
                                 'username': 'username',
                                 'password': 'pwd',
                                 'port': 'port'},
                                'log':
                                {'file': 'data/taca.log'}}
        self.assertItemsEqual(expexted_config_data, got_config_data)
        with self.assertRaises(IOError):
            missing_config_data = config.load_yaml_config('data/missing_file.yaml)')


class TestBioinfoTab(unittest.TestCase):
    """ Test class for bioinfo_tab """

    @classmethod
    def setUpClass(self):
        self.rootdir = tempfile.mkdtemp(prefix="test_taca_bt")
        self.new_run = os.path.join(self.rootdir,'nosync/190821_M01545_0252_000000001')
        os.makedirs(self.new_run)
        self.demux_run = os.path.join(self.rootdir, '190821_M01545_0252_000000002')
        os.makedirs(os.path.join(self.demux_run, 'Unaligned_1'))
        self.seq_run = os.path.join(self.rootdir, '190821_M01545_0252_000000003')
        os.makedirs(self.seq_run)
        self.error_run = os.path.join(self.rootdir, '190821_M01545_0252_000000004')
        os.makedirs(self.error_run)
        with open(os.path.join(self.error_run, 'RTAComplete.txt'), 'w') as fh:
            fh.write("This is some contents\n")

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.rootdir)

    def test_setupServer(self):
        """ Set up server connection """
        config = {'statusdb':
                  {'url': 'url',
                   'username': 'username',
                   'password': 'pwd',
                   'port': '1234'}}
        got_connection = bioinfo_tab.setupServer(config)
        assert isinstance(got_connection, couchdb.Server)

    @mock.patch('taca.utils.bioinfo_tab.update_statusdb', return_value=None)
    def test_collect_runs(self, mock_update_statusdb):
        """ find runs in specified dir """
        bioinfo_tab.collect_runs()
        calls = [mock.call('data/test_data/190201_A00621_0032_BHHFCFDSXX'), mock.call('data/test_data/nosync/190201_A00621_0032_BHHFCFDSXY')]
        mock_update_statusdb.assert_has_calls(calls)

    def test_get_status_new(self):
        """ return status New """
        got_status = bioinfo_tab.get_status(self.new_run)
        self.assertEqual(got_status, 'New')

    def test_get_status_demultiplexing(self):
        """ return status Demultiplexing """
        got_status = bioinfo_tab.get_status(self.demux_run)
        self.assertEqual(got_status, 'Demultiplexing')

    def test_get_status_sequencing(self):
        """ return status Sequencing """
        got_status = bioinfo_tab.get_status(self.seq_run)
        self.assertEqual(got_status, 'Sequencing')

    def test_get_status_error(self):
        """ return status ERROR """
        got_status = bioinfo_tab.get_status(self.error_run)
        self.assertEqual(got_status, 'ERROR')

    def test_parse_sample_sheet(self):
        """ parse samplesheet """
        sample_sheet = 'data/samplesheet.csv'
        expected_data = [{'SampleWell': '1:1',
                          'index': 'GAATTCGT',
                          'Lane': '1',
                          'SamplePlate': 'FCB_150423',
                          'SampleName': 'P1775_147',
                          'SampleID': 'Sample_P1775_147',
                          'Project': 'J_Lundeberg_14_24'}]
        parsed_data = bioinfo_tab.parse_samplesheet(sample_sheet, "run_dir")
        self.assertEqual(expected_data, parsed_data)

    def test_parse_sample_sheet_is_miseq(self):
        """ parse miseq samplesheet """
        sample_sheet = 'data/miseq_samplesheet.csv'
        expected_data = [{'SampleWell': '1:1',
                          'index': 'GAATTCGT',
                          'Lane': '1',
                          'SamplePlate': 'FCB_150423',
                          'SampleName': 'P1775_147',
                          'SampleID': 'Sample_P1775_147',
                          'Project': 'J_Lundeberg_14_24'}]
        parsed_data = bioinfo_tab.parse_samplesheet(sample_sheet, "run_dir", is_miseq=True)
        self.assertEqual(expected_data, parsed_data)

    def test_parse_sample_sheet_is_miseq_error(self):
        """ return empty list if not production or application in miseq sample sheet """
        sample_sheet = 'data/samplesheet.csv'
        parsed_data = bioinfo_tab.parse_samplesheet(sample_sheet, "run_dir", is_miseq=True)
        self.assertEqual(parsed_data, [])

    @mock.patch('taca.utils.bioinfo_tab.send_mail')
    @mock.patch('taca.utils.bioinfo_tab.datetime.datetime')
    def test_error_mailer_no_samplesheet(self, mock_datetime, mock_send_mail):
        """ send email if no_samplesheet error """
        body='TACA has encountered an issue that might be worth investigating\n'
        body+='The offending entry is: '
        body+= 'run_missing_samplesheet'
        body+='\n\nSincerely, TACA'
        subject='ERROR, Samplesheet error'
        mock_datetime.now()
        mock_datetime.now().hour = 7
        bioinfo_tab.error_emailer('no_samplesheet', 'run_missing_samplesheet')
        mock_send_mail.assert_called_with(subject, body, 'some_user@some_email.com')

    @mock.patch('taca.utils.bioinfo_tab.send_mail')
    @mock.patch('taca.utils.bioinfo_tab.datetime.datetime')
    def test_error_mailer_failed_run(self, mock_datetime, mock_send_mail):
        """ send email if failed_run error """
        body='TACA has encountered an issue that might be worth investigating\n'
        body+='The offending entry is: '
        body+= 'failed_run'
        body+='\n\nSincerely, TACA'
        subject='WARNING, Reinitialization of partially failed FC'
        mock_datetime.now()
        mock_datetime.now().hour = 7
        bioinfo_tab.error_emailer('failed_run', 'failed_run')
        mock_send_mail.assert_called_with(subject, body, 'some_user@some_email.com')

    @mock.patch('taca.utils.bioinfo_tab.send_mail')
    @mock.patch('taca.utils.bioinfo_tab.datetime.datetime')
    def test_error_mailer_weird_samplesheet(self, mock_datetime, mock_send_mail):
        """ send email if weird_samplesheet error """
        body='TACA has encountered an issue that might be worth investigating\n'
        body+='The offending entry is: '
        body+= 'weird_samplesheet_run'
        body+='\n\nSincerely, TACA'
        subject='ERROR, Incorrectly formatted samplesheet'
        mock_datetime.now()
        mock_datetime.now().hour = 7
        bioinfo_tab.error_emailer('weird_samplesheet', 'weird_samplesheet_run')
        mock_send_mail.assert_called_with(subject, body, 'some_user@some_email.com')
