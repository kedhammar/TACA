"""Backup methods and utilities"""
import couchdb
import logging
import os
import re
import sys
import shutil
import subprocess as sp

from datetime import datetime
from taca.utils.config import CONFIG
from taca.utils import filesystem, misc

logger = logging.getLogger(__name__)

class run_vars(object):
    """A simple variable storage class"""
    def __init__(self, run):
        self.abs_path = os.path.abspath(run)
        self.path, self.name = os.path.split(self.abs_path)
        self.name = self.name.split('.', 1)[0]
        self.zip = "{}.tar.gz".format(self.name)
        self.key = "{}.key".format(self.name)
        self.key_encrypted = "{}.key.gpg".format(self.name)
        self.zip_encrypted = "{}.tar.gz.gpg".format(self.name)

class backup_utils(object):
    """A class object with main utility methods related to backing up"""

    def __init__(self, run=None):
        self.run = run
        self.fetch_config_info()

    def fetch_config_info(self):
        """Try to fecth required info from the config file. Log and exit if any neccesary info is missing"""
        try:
            self.data_dirs = CONFIG['backup']['data_dirs']
            self.archive_dirs = CONFIG['backup']['archive_dirs']
            self.keys_path = CONFIG['backup']['keys_path']
            self.gpg_receiver = CONFIG['backup']['gpg_receiver']
            self.couch_info = CONFIG['statusdb']
        except KeyError as e:
            logger.error("Config file is missing the key {}, make sure it have all required information".format(str(e)))
            raise SystemExit

    def collect_runs(self, ext=None, filter_by_ext=False):
        """Collect runs from archive directeries"""
        self.runs = []
        if self.run:
            run = run_vars(self.run)
            if not re.match(filesystem.RUN_RE, run.name):
                logger.error("Given run {} did not match a FC pattern".format(self.run))
                raise SystemExit
            self.runs.append(run)
        else:
            for adir in self.archive_dirs.values():
                if not os.path.isdir(adir):
                    logger.warn("Path {} does not exist or it is not a directory".format(adir))
                    return self.runs
                for item in os.listdir(adir):
                    if filter_by_ext and not item.endswith(ext):
                        continue
                    elif item.endswith(ext):
                        item = item.replace(ext, '')
                    elif not os.path.isdir(item):
                        continue
                    if re.match(filesystem.RUN_RE, item) and item not in self.runs:
                        self.runs.append(run_vars(os.path.join(adir, item)))

    def run_is_demuxed(self, run):
        """Check in StatusDB 'x_flowcells' database if the given run has an entry which means it was
        demultiplexed (as TACA only creates a document upon successfull demultiplexing)"""
        run_terms = run.split('_')
        run_date = run_terms[0]
        run_fc = run_terms[-1]
        run_name = "{}_{}".format(run_date, run_fc)
        couch_info = self.couch_info
        # connect to statusdb using info fectched from config file
        try:
            server = "http://{username}:{password}@{url}:{port}".format(url=couch_info['url'],username=couch_info['username'],
                                                                        password=couch_info['password'],port=couch_info['port'])
            couch = couchdb.Server(server)
            fc_db = couch[couch_info['db']]
            fc_names = [entry.key for entry in fc_db.view("names/name", reduce=False)]
        except Exception, e:
            logger.error(e.message)
            raise
        if run_name in fc_names:
            logger.info("Run {} is demultiplexed and proceeding with encryption".format(run))
            return True
        else:
            logger.warn("Run {} is not demultiplexed yet, so skipping it".format(run))
            return False

    def avail_disk_space(self, path, run):
        """Check the space on file system based on parent directory of the run"""
        # not able to fetch runtype use the max size as precaution, size units in GB
        illumina_run_sizes = {'hiseq' : 500, 'hiseqx' : 900, 'miseq' : 20}
        required_size = illumina_run_sizes.get(self._get_run_type(run), 900) * 2
        # check for any ongoing runs and add up the required size accrdingly
        for ddir in self.data_dirs.values():
            for item in os.listdir(ddir):
                if not re.match(filesystem.RUN_RE, item):
                    continue
                if not os.path.exists(os.path.join(ddir, item, "RTAComplete.txt")):
                    required_size += illumina_run_sizes.get(self._get_run_type(run), 900)
        # get available free space from the file system
        try:
            df_proc = sp.Popen(['df', path], stdout=sp.PIPE, stderr=sp.PIPE)
            df_out, df_err = df_proc.communicate()
            available_size = int(df_out.strip().split('\n')[-1].strip().split()[2])/1024/1024
        except Exception, e:
            logger.error("Evaluation of disk space failed with error {}".format(e))
            return False
        if available_size > required_size:
            return True
        return False

    def _get_run_type(self, run):
        """Returns run type based on the flowcell name"""
        run_type = ''
        try:
            if "ST-" in run:
                run_type = "hiseqx"
            elif "-" in run.split('_')[-1]:
                run_type = "miseq"
            else:
                run_type = "hiseq"
        except:
            logger.warn("Could not fetch run type for run {}".format(run))
        return run_type

    def _call_commands(self, cmd1, cmd2=None, out_file=None, return_out=False, tmp_files=[]):
        """Call an external command(s) with atmost two commands per function call.
        Given 'out_file' is always used for the later cmd and also stdout can be return
        for the later cmd. In case of failure, the 'tmp_files' are removed"""
        if out_file:
            if not cmd2:
                stdout1 = open(out_file, 'w')
            else:
                stdout1 = sp.PIPE
                stdout2 = open(out_file, 'w')
        else:
            stdout1 = sp.PIPE
            stdout2 = sp.PIPE
        # calling the commands
        try:
            cmd1 = cmd1.split()
            p1 = sp.Popen(cmd1, stdout=stdout1, stderr=sp.PIPE)
            if cmd2:
                cmd2 = cmd2.split()
                p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=stdout2, stderr=sp.PIPE)
                p2_stat = p2.wait()
                p2_out, p2_err = p2.communicate()
                if not self._check_status(cmd2, p2_stat, p2_err, tmp_files):
                    return False
            p1_stat = p1.wait()
            p1_out, p1_err = p1.communicate()
            if not self._check_status(cmd1, p1_stat, p1_err, tmp_files):
                return False
        except Exception, e:
            raise e
        finally:
            if out_file:
                if not cmd2:
                    stdout1.close()
                else:
                    stdout2.close()
            if return_out:
                if cmd2:
                    return (True, p2_out)
                return (True, p1_out)
            return True

    def _check_status(self, cmd, status, err_msg, files_to_remove=[]):
        """Check if a subprocess status is success and log error if failed"""
        if status != 0:
            self._clean_tmp_files(files_to_remove)
            logger.error("Command '{}' failed with the error '{}'".format(" ".join(cmd),err_msg))
            return False
        return True

    def _clean_tmp_files(self, files):
        """Remove the file is exist"""
        for fl in files:
            if os.path.exists(fl):
                os.remove(fl)

    @classmethod
    def encrypt_runs(cls, run, force):
        """Encrypt the runs that have been collected"""
        bk = cls(run)
        bk.collect_runs(ext=".tar.gz")
        logger.info("In total, found {} run(s) to be encrypted".format(len(bk.runs)))
        for run in bk.runs:
            run.flag = "{}.encrypting".format(run.name)
            run.dst_key_encrypted = os.path.join(bk.keys_path, run.key_encrypted)
            tmp_files = [run.zip_encrypted, run.key_encrypted, run.key, run.flag]
            logger.info("Encryption of run {} is now started".format(run.name))
            # Check if there is enough space
            if not bk.avail_disk_space(path=run.path, run=run.name):
                logger.error("There is no enough disk space for compression, kindly check and archive encrypted runs")
                raise SystemExit
            # Check if the run in demultiplexed
            if not force and not bk.run_is_demuxed(run.name):
                continue
            with filesystem.chdir(run.path):
                # skip run if already ongoing
                if os.path.exists(run.flag):
                    logger.warn("Run {} is already being encrypted, so skipping now".format(run.name))
                    continue
                flag = open(run.flag, 'w').close()
                # zip the run directory
                if os.path.exists(run.zip):
                    if os.path.isdir(run.name):
                        logger.warn("Both run source and zipped archive exist for run {}, skipping run as precaution".format(run.name))
                        bk._clean_tmp_files([run.flag])
                        continue
                    logger.info("Zipped archive already exist for run {}, so using it for encryption".format(run.name))
                else:
                    logger.info("Creating zipped archive for run {}".format(run.name))
                    if bk._call_commands(cmd1="tar -cf - {}".format(run.name), cmd2="pigz --fast -c -",
                                         out_file=run.zip, tmp_files=[run.zip, run.flag]):
                        logger.info("Run {} was successfully compressed, so removing the run source directory".format(run.name))
                        shutil.rmtree(run.name)
                    else:
                        logger.warn("Skipping run {} and moving on".format(run.name))
                        continue
                # Remove encrypted file if already exists
                if os.path.exists(run.zip_encrypted):
                    logger.warn(("Removing already existing encrypted file for run {}, this is a precaution "
                                 "to make sure the file was encrypted with correct key file".format(run.name)))
                    bk._clean_tmp_files([run.zip_encrypted, run.key, run.key_encrypted, run.dst_key_encrypted])
                # Generate random key to use as pasphrase
                if not bk._call_commands(cmd1="gpg --gen-random 1 256", out_file=run.key, tmp_files=tmp_files):
                    logger.warn("Skipping run {} and moving on".format(run.name))
                    continue
                logger.info("Generated randon phrase key for run {}".format(run.name))
                # Calculate md5 sum pre encryption
                if not force:
                    logger.info("Calculating md5sum before encryption")
                    md5_call, md5_out = bk._call_commands(cmd1="md5sum {}".format(run.zip), return_out=True, tmp_files=tmp_files)
                    if not md5_call:
                        logger.warn("Skipping run {} and moving on".format(run.name))
                        continue
                    md5_pre_encrypt = md5_out.split()[0]
                # Encrypt the zipped run file
                logger.info("Encrypting the zipped run file")
                if not bk._call_commands(cmd1=("gpg --symmetric --cipher-algo aes256 --passphrase-file {} --batch --compress-algo "
                                               "none -o {} {}".format(run.key, run.zip_encrypted, run.zip)), tmp_files=tmp_files):
                    logger.warn("Skipping run {} and moving on".format(run.name))
                    continue
                # Decrypt and check for md5
                if not force:
                    logger.info("Calculating md5sum after encryption")
                    md5_call, md5_out = bk._call_commands(cmd1="gpg --decrypt --cipher-algo aes256 --passphrase-file {} --batch {}".format(run.key, run.zip_encrypted),
                                                          cmd2="md5sum", return_out=True, tmp_files=tmp_files)
                    if not md5_call:
                        logger.warn("Skipping run {} and moving on".format(run.name))
                        continue
                    md5_post_encrypt = md5_out.split()[0]
                    if md5_pre_encrypt != md5_post_encrypt:
                        logger.error(("md5sum did not match before {} and after {} encryption. Will remove temp files and "
                                      "move on".format(md5_pre_encrypt, md5_post_encrypt)))
                        bk._clean_tmp_files(tmp_files)
                        continue
                    logger.info("Md5sum is macthing before and after encryption")
                # Encrypt and move the key file
                if bk._call_commands(cmd1="gpg -e -r {} -o {} {}".format(bk.gpg_receiver, run.key_encrypted, run.key), tmp_files=tmp_files):
                    shutil.move(run.key_encrypted, run.dst_key_encrypted)
                else:
                    logger.error("Encrption of key file failed, skipping run")
                    continue
                bk._clean_tmp_files([run.zip, run.key, run.flag])
                logger.info("Encryption of run {} is successfully done, removing zipped run file".format(run.name))

    @classmethod
    def pdc_put(cls, run):
        """Archive the collected runs to PDC"""
        bk = cls(run)
        bk.collect_runs(ext=".tar.gz.gpg", filter_by_ext=True)
        logger.info("In total, found {} run(s) to send PDC".format(len(bk.runs)))
        for run in bk.runs:
            run.flag = "{}.archiving".format(run.name)
            run.dst_key_encrypted = os.path.join(bk.keys_path, run.key_encrypted)
            if run.path not in bk.archive_dirs.values():
                logger.error(("Given run is not in one of the archive directories {}. Kindly move the run {} to appropriate "
                              "archive dir before sending it to PDC".format(",".join(bk.archive_dirs.values()), run.name)))
                continue
            if not os.path.exists(run.dst_key_encrypted):
                logger.error("Encrypted key file {} is not found for file {}, skipping it".format(run.dst_key_encrypted, run.zip_encrypted))
                continue
            # skip run if already ongoing
            if os.path.exists(run.flag):
                logger.warn("Run {} is already being archived, so skipping now".format(run.name))
                continue
            flag = open(run.flag, 'w').close()
            with filesystem.chdir(run.path):
                if bk._call_commands(cmd1="dsmc archive {}".format(run.zip_encrypted), tmp_files=[run.flag]):
                    if bk._call_commands(cmd1="dsmc archive {}".format(run.dst_key_encrypted), tmp_files=[run.flag]):
                        logger.info("Successfully sent file {} to PDC, removing file locally from {}".format(run.zip_encrypted, run.path))
                        shutil.move(run.zip_encrypted, os.path.join("sent_data", run.zip_encrypted))
                        shutil.move(run.dst_key_encrypted, os.path.join("sent_data", run.key_encrypted))
                        continue
                logger.warn("Sending file {} to PDC failed".format(run.zip_encrypted))

