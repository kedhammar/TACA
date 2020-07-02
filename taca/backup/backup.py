"""Backup methods and utilities."""
import couchdb
import logging
import os
import re
import shutil
import subprocess as sp
import time

from datetime import datetime
from taca.utils.config import CONFIG
from taca.utils import filesystem, misc

logger = logging.getLogger(__name__)

class run_vars(object):
    """A simple variable storage class."""
    def __init__(self, run):
        self.abs_path = os.path.abspath(run)
        self.path, self.name = os.path.split(self.abs_path)
        self.name = self.name.split('.', 1)[0]
        self.zip = '{}.tar.gz'.format(self.name)
        self.key = '{}.key'.format(self.name)
        self.key_encrypted = '{}.key.gpg'.format(self.name)
        self.zip_encrypted = '{}.tar.gz.gpg'.format(self.name)

class backup_utils(object):
    """A class object with main utility methods related to backing up."""

    def __init__(self, run=None):
        self.run = run
        self.fetch_config_info()
        self.host_name = os.getenv('HOSTNAME', os.uname()[1]).split('.', 1)[0]

    def fetch_config_info(self):
        """Try to fecth required info from the config file. Log and exit if any neccesary info is missing."""
        try:
            self.data_dirs = CONFIG['backup']['data_dirs']
            self.archive_dirs = CONFIG['backup']['archive_dirs']
            self.keys_path = CONFIG['backup']['keys_path']
            self.gpg_receiver = CONFIG['backup']['gpg_receiver']
            self.mail_recipients = CONFIG['mail']['recipients']
            self.check_demux = CONFIG.get('backup', {}).get('check_demux', False)
            self.couch_info = CONFIG.get('statusdb')
        except KeyError as e:
            logger.error('Config file is missing the key {}, make sure it have all required information'.format(str(e)))
            raise SystemExit

    def collect_runs(self, ext=None, filter_by_ext=False):
        """Collect runs from archive directories."""
        self.runs = []
        if self.run:
            run = run_vars(self.run)
            if not re.match(filesystem.RUN_RE, run.name):
                logger.error('Given run {} did not match a FC pattern'.format(self.run))
                raise SystemExit
            self.runs.append(run)
        else:
            for adir in self.archive_dirs.values():
                if not os.path.isdir(adir):
                    logger.warn('Path {} does not exist or it is not a directory'.format(adir))
                    continue
                for item in os.listdir(adir):
                    if filter_by_ext and not item.endswith(ext):
                        continue
                    elif item.endswith(ext):
                        item = item.replace(ext, '')
                    elif not os.path.isdir(os.path.join(adir, item)):
                        continue
                    if re.match(filesystem.RUN_RE, item) and item not in self.runs:
                        self.runs.append(run_vars(os.path.join(adir, item)))

    def avail_disk_space(self, path, run):
        """Check the space on file system based on parent directory of the run."""
        # not able to fetch runtype use the max size as precaution, size units in GB
        illumina_run_sizes = {'hiseq': 500, 'hiseqx': 900, 'novaseq': 1800, 'miseq': 20, 'nextseq': 50}
        required_size = illumina_run_sizes.get(self._get_run_type(run), 900) * 2
        # check for any ongoing runs and add up the required size accrdingly
        for ddir in self.data_dirs.values():
            if not os.path.isdir(ddir):
                continue
            for item in os.listdir(ddir):
                if not re.match(filesystem.RUN_RE, item):
                    continue
                if not os.path.exists(os.path.join(ddir, item, 'RTAComplete.txt')):
                    required_size += illumina_run_sizes.get(self._get_run_type(run), 900)
        # get available free space from the file system
        try:
            df_proc = sp.Popen(['df', path], stdout=sp.PIPE, stderr=sp.PIPE)
            df_out, df_err = df_proc.communicate()
            available_size = int(df_out.strip().split('\n')[-1].strip().split()[2])/1024/1024
        except Exception as e:
            logger.error('Evaluation of disk space failed with error {}'.format(e))
            raise SystemExit
        if available_size < required_size:
            e_msg = 'Required space for encryption is {}GB, but only {}GB available'.format(required_size, available_size)
            subjt = 'Low space for encryption - {}'.format(self.host_name)
            logger.error(e_msg)
            misc.send_mail(subjt, e_msg, self.mail_recipients)
            raise SystemExit

    def file_in_pdc(self, src_file, silent=True):
        """Check if the given files exist in PDC."""
        # dsmc will return zero/True only when file exists, it returns
        # non-zero/False though cmd is execudted but file not found
        src_file_abs = os.path.abspath(src_file)
        try:
            sp.check_call(['dsmc', 'query', 'archive', src_file_abs], stdout=sp.PIPE, stderr=sp.PIPE)
            value = True
        except sp.CalledProcessError:
            value = False
        if not silent:
            msg = 'File {} {} in PDC'.format(src_file_abs, 'exist' if value else 'do not exist')
            logger.info(msg)
        return value

    def _get_run_type(self, run):
        """Returns run type based on the flowcell name."""
        run_type = ''
        try:
            if 'ST-' in run:
                run_type = 'hiseqx'
            elif '-' in run.split('_')[-1]:
                run_type = 'miseq'
            elif '_A0' in run:
                run_type = 'novaseq'
            elif '_NS' in run:
                run_type = 'nextseq'
            else:
                run_type = 'hiseq'
        except:
            logger.warn('Could not fetch run type for run {}'.format(run))
        return run_type

    def _call_commands(self, cmd1, cmd2=None, out_file=None, return_out=False, mail_failed=False, tmp_files=[]):
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
                if not self._check_status(cmd2, p2_stat, p2_err, mail_failed, tmp_files):
                    return (False, p2_err) if return_out else False
            p1_stat = p1.wait()
            p1_out, p1_err = p1.communicate()
            if not self._check_status(cmd1, p1_stat, p1_err, mail_failed, tmp_files):
                return (False, p1_err) if return_out else False
            if return_out:
                return (True, p2_out) if cmd2 else (True, p1_out)
            return True
        except Exception as e:
            raise e
        finally:
            if out_file:
                if not cmd2:
                    stdout1.close()
                else:
                    stdout2.close()

    def _check_status(self, cmd, status, err_msg, mail_failed, files_to_remove=[]):
        """Check if a subprocess status is success and log error if failed."""
        if status != 0:
            self._clean_tmp_files(files_to_remove)
            if mail_failed:
                subjt = 'Command call failed - {}'.format(self.host_name)
                e_msg = 'Called cmd: {}\n\nError msg: {}'.format(' '.join(cmd), err_msg)
                misc.send_mail(subjt, e_msg, self.mail_recipients)
            logger.error('Command "{}" failed with the error "{}"'.format(' '.join(cmd),err_msg))
            return False
        return True

    def _clean_tmp_files(self, files):
        """Remove the file is exist."""
        for fl in files:
            if os.path.exists(fl):
                os.remove(fl)

    def _log_pdc_statusdb(self, run):
        """Log the time stamp in statusDB if a file is succussfully sent to PDC."""
        try:
            run_vals = run.split('_')
            run_fc = '{}_{}'.format(run_vals[0],run_vals[-1])
            server = 'http://{username}:{password}@{url}:{port}'.format(url=self.couch_info['url'],username=self.couch_info['username'],
                                                                        password=self.couch_info['password'],port=self.couch_info['port'])
            couch = couchdb.Server(server)
            db = couch[self.couch_info['db']]
            fc_names = {e.key:e.id for e in db.view('names/name', reduce=False)}
            d_id = fc_names[run_fc]
            doc = db.get(d_id)
            doc['pdc_archived'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            db.save(doc)
            logger.info('Logged "pdc_archived" timestamp for fc {} in statusdb doc "{}"'.format(run, d_id))
        except:
            logger.warn('Not able to log "pdc_archived" timestamp for run {}'.format(run))

    @classmethod
    def encrypt_runs(cls, run, force):
        """Encrypt the runs that have been collected."""
        bk = cls(run)
        bk.collect_runs(ext='.tar.gz')
        logger.info('In total, found {} run(s) to be encrypted'.format(len(bk.runs)))
        for run in bk.runs:
            run.flag = '{}.encrypting'.format(run.name)
            run.dst_key_encrypted = os.path.join(bk.keys_path, run.key_encrypted)
            tmp_files = [run.zip_encrypted, run.key_encrypted, run.key, run.flag]
            logger.info('Encryption of run {} is now started'.format(run.name))
            # Check if there is enough space and exit if not
            bk.avail_disk_space(run.path, run.name)
            # Check if the run in demultiplexed
            if not force and bk.check_demux:
                if not misc.run_is_demuxed(run.name, bk.couch_info):
                    logger.warn('Run {} is not demultiplexed yet, so skipping it'.format(run.name))
                    continue
                logger.info('Run {} is demultiplexed and proceeding with encryption'.format(run.name))
            with filesystem.chdir(run.path):
                # skip run if already ongoing
                if os.path.exists(run.flag):
                    logger.warn('Run {} is already being encrypted, so skipping now'.format(run.name))
                    continue
                flag = open(run.flag, 'w').close()
                # zip the run directory
                if os.path.exists(run.zip):
                    if os.path.isdir(run.name):
                        logger.warn('Both run source and zipped archive exist for run {}, skipping run as precaution'.format(run.name))
                        bk._clean_tmp_files([run.flag])
                        continue
                    logger.info('Zipped archive already exist for run {}, so using it for encryption'.format(run.name))
                else:
                    logger.info('Creating zipped archive for run {}'.format(run.name))
                    if bk._call_commands(cmd1='tar -cf - {}'.format(run.name), cmd2='pigz --fast -c -',
                                         out_file=run.zip, mail_failed=True, tmp_files=[run.zip, run.flag]):
                        logger.info('Run {} was successfully compressed, so removing the run source directory'.format(run.name))
                        shutil.rmtree(run.name)
                    else:
                        logger.warn('Skipping run {} and moving on'.format(run.name))
                        continue
                # Remove encrypted file if already exists
                if os.path.exists(run.zip_encrypted):
                    logger.warn(('Removing already existing encrypted file for run {}, this is a precaution '
                                 'to make sure the file was encrypted with correct key file'.format(run.name)))
                    bk._clean_tmp_files([run.zip_encrypted, run.key, run.key_encrypted, run.dst_key_encrypted])
                # Generate random key to use as pasphrase
                if not bk._call_commands(cmd1='gpg --gen-random 1 256', out_file=run.key, tmp_files=tmp_files):
                    logger.warn('Skipping run {} and moving on'.format(run.name))
                    continue
                logger.info('Generated randon phrase key for run {}'.format(run.name))
                # Calculate md5 sum pre encryption
                if not force:
                    logger.info('Calculating md5sum before encryption')
                    md5_call, md5_out = bk._call_commands(cmd1='md5sum {}'.format(run.zip), return_out=True, tmp_files=tmp_files)
                    if not md5_call:
                        logger.warn('Skipping run {} and moving on'.format(run.name))
                        continue
                    md5_pre_encrypt = md5_out.split()[0]
                # Encrypt the zipped run file
                logger.info('Encrypting the zipped run file')
                if not bk._call_commands(cmd1=('gpg --symmetric --cipher-algo aes256 --passphrase-file {} --batch --compress-algo '
                                               'none -o {} {}'.format(run.key, run.zip_encrypted, run.zip)), tmp_files=tmp_files):
                    logger.warn('Skipping run {} and moving on'.format(run.name))
                    continue
                # Decrypt and check for md5
                if not force:
                    logger.info('Calculating md5sum after encryption')
                    md5_call, md5_out = bk._call_commands(cmd1='gpg --decrypt --cipher-algo aes256 --passphrase-file {} --batch {}'.format(run.key, run.zip_encrypted),
                                                          cmd2='md5sum', return_out=True, tmp_files=tmp_files)
                    if not md5_call:
                        logger.warn('Skipping run {} and moving on'.format(run.name))
                        continue
                    md5_post_encrypt = md5_out.split()[0]
                    if md5_pre_encrypt != md5_post_encrypt:
                        logger.error(('md5sum did not match before {} and after {} encryption. Will remove temp files and '
                                      'move on'.format(md5_pre_encrypt, md5_post_encrypt)))
                        bk._clean_tmp_files(tmp_files)
                        continue
                    logger.info('Md5sum is macthing before and after encryption')
                # Encrypt and move the key file
                if bk._call_commands(cmd1='gpg -e -r {} -o {} {}'.format(bk.gpg_receiver, run.key_encrypted, run.key), tmp_files=tmp_files):
                    shutil.move(run.key_encrypted, run.dst_key_encrypted)
                else:
                    logger.error('Encrption of key file failed, skipping run')
                    continue
                bk._clean_tmp_files([run.zip, run.key, run.flag])
                logger.info('Encryption of run {} is successfully done, removing zipped run file'.format(run.name))

    @classmethod
    def pdc_put(cls, run):
        """Archive the collected runs to PDC."""
        bk = cls(run)
        bk.collect_runs(ext='.tar.gz.gpg', filter_by_ext=True)
        logger.info('In total, found {} run(s) to send PDC'.format(len(bk.runs)))
        for run in bk.runs:
            run.flag = '{}.archiving'.format(run.name)
            run.dst_key_encrypted = os.path.join(bk.keys_path, run.key_encrypted)
            if run.path not in bk.archive_dirs.values():
                logger.error(('Given run is not in one of the archive directories {}. Kindly move the run {} to appropriate '
                              'archive dir before sending it to PDC'.format(','.join(bk.archive_dirs.values()), run.name)))
                continue
            if not os.path.exists(run.dst_key_encrypted):
                logger.error('Encrypted key file {} is not found for file {}, skipping it'.format(run.dst_key_encrypted, run.zip_encrypted))
                continue
            with filesystem.chdir(run.path):
                #skip run if being encrypted
                if os.path.exists('{}.encrypting'.format(run.name)):
                    logger.warn('Run {} is currently being encrypted, so skipping now'.format(run.name))
                    continue
                # skip run if already ongoing
                if os.path.exists(run.flag):
                    logger.warn('Run {} is already being archived, so skipping now'.format(run.name))
                    continue
                if bk.file_in_pdc(run.zip_encrypted, silent=False) or bk.file_in_pdc(run.dst_key_encrypted, silent=False):
                    logger.warn('Seems like files realted to run {} already exist in PDC, check and cleanup'.format(run.name))
                    continue
                flag = open(run.flag, 'w').close()
                logger.info('Sending file {} to PDC'.format(run.zip_encrypted))
                if bk._call_commands(cmd1='dsmc archive {}'.format(run.zip_encrypted), tmp_files=[run.flag]):
                    time.sleep(15) # give some time just in case 'dsmc' needs to settle
                    if bk._call_commands(cmd1='dsmc archive {}'.format(run.dst_key_encrypted), tmp_files=[run.flag]):
                        time.sleep(5) # give some time just in case 'dsmc' needs to settle
                        if bk.file_in_pdc(run.zip_encrypted) and bk.file_in_pdc(run.dst_key_encrypted):
                            logger.info('Successfully sent file {} to PDC, removing file locally from {}'.format(run.zip_encrypted, run.path))
                            if bk.couch_info:
                                bk._log_pdc_statusdb(run.name)
                            bk._clean_tmp_files([run.zip_encrypted, run.dst_key_encrypted, run.flag])
                        continue
                logger.warn('Sending file {} to PDC failed'.format(run.zip_encrypted))
