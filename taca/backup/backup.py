"""Backup methods and utilities"""
import couchdb
import logging
import os
import re
import sys
import shutil
import subprocess as sp

from datetime import datetime
from time import sleep
from taca.utils.config import CONFIG
from taca.utils import filesystem, misc, statusdb

logger = logging.getLogger(__name__)

def encrypt_data(run, force):
    """Compress and encrypt the givne run if given or collect runs from
    specifed root in config file
    
    :param str run: optional to give particular run for encryption
    :param boolean force: boolean to say to ignore basic checks
    """
    #try fetchig required info from CONFIG file, complain if couldn't
    try:
        data_dirs = CONFIG['backup']['data_dirs']
        archive_dirs = CONFIG['backup']['archive_dirs']
        keys_path = CONFIG['backup']['keys_path']
        couch_info = CONFIG['statusdb']
    except KeyError as e:
        logger.error("Config file is missing the key {}, make sure it have all required information".format(str(e)))
        raise SystemExit
    
    runs = []
    if run:
        run_tmp = run.split('.', 1)[0] # in case if the given run was a zipped file
        if not re.match(filesystem.RUN_RE, os.path.basename(run_tmp)):
            logger.error("Given run {} did not match a FC pattern".format(run))
            raise SystemExit
        runs.append(run_tmp)
    else:
        for adir in archive_dirs:
            runs.extend(_get_runs_from_path(path=adir, ext=".tar.gz"))
    
    logger.info("In total, found {} run(s) to be encrypted".format(len(runs)))
    for run in runs:
        run_path = os.path.dirname(run)
        run = os.path.basename(run)
        run_zip = "{}.tar.gz".format(run)
        run_key = "{}.key".format(run)
        run_flag = "{}.encrypting".format(run)
        run_key_encrypted = "{}.key.gpg".format(run)
        run_zip_encrypted = "{}.tar.gz.gpg".format(run)
        dst_key_encrypted = os.path.join(keys_path, run_key_encrypted)
        _tmp_files = [run_zip_encrypted, run_key_encrypted, run_key, run_flag]
        logger.info("Encryption of run {} is now started".format(run))
        # Check if there is enough space 
        if not force and not _avail_disk_space(path=run_path, run=run, data_dirs=data_dirs):
            logger.error("There is no enough disk space for compression, kindly check and archive encrypted runs")
            raise SystemExit
        # Check if the run in demultiplexed
        if not force and not _run_is_demuxed(run, couch_info):
            continue
        
        with filesystem.chdir(run_path):
            # skip run if already ongoing
            if os.path.exists(run_flag):
                logger.warn("Run {} is already being encrypted, so skipping now".format(run))
                continue
            flag = open(run_flag, 'w').close()
            
            # zip the run directory
            if os.path.exists(run_zip):
                if os.path.isdir(run):
                    logger.warn("Both run source and zipped archive exist for run {}, skipping run as precaution".format(run))
                    _clean_tmp_files([run_flag])
                    continue
                logger.info("Zipped archive already exist for run {}, so using it for encryption".format(run))
            else:
                logger.info("Creating zipped archive for run {}".format(run))
                with open(run_zip, 'w') as rzip:
                    tar_proc = sp.Popen(['tar', '-cf', '-', run], stdout=sp.PIPE, stderr=sp.PIPE)
                    pigz_proc = sp.Popen(['pigz', '--fast', '-c', '-'], stdin=tar_proc.stdout, stdout=rzip, stderr=sp.PIPE)
                    tar_status = tar_proc.wait()
                    pigz_status = pigz_proc.wait()
                # Check if the process have completed succesfully, Popen returns 0 for successfull completion
                if tar_status == 0 and pigz_status == 0:
                    logger.info("Run {} was successfully compressed, so removing the run source directory".format(run))
                    shutil.rmtree(run)
                else:
                    if os.path.exists(run_zip):
                        os.remove(run_zip)
                    err_msg = "Compresion for run {} failed with errors, so skipping it now.".format(run)
                    if tar_status:
                        err_msg = "{}\n{}".format(err_msg, tar_proc.communicate()[-1])
                    if pigz_status:
                        err_msg = "{}\n{}".format(err_msg, pigz_proc.communicate()[-1])
                    logger.error(err_msg)
                    _clean_tmp_files(_tmp_files)
                    continue

            # Remove encrypted file if already exists
            if os.path.exists(run_zip_encrypted):
                logger.warn(("Removing already existing encrypted file for run {}, this is a precaution "
                             "to make sure the file was encrypted with correct key file".format(run)))
                _clean_tmp_files([run_zip_encrypted, run_key, run_key_encrypted, dst_key_encrypted])
            
            # Generate random key to use as pasphrase:
            with open(run_key, 'w') as kfl:
                rkey_proc = sp.Popen(['gpg', '--gen-random', '1', '256'], stdout=kfl, stderr=sp.PIPE)
                rkey_status = rkey_proc.wait()
                rkey_out, rkey_err = rkey_proc.communicate()
            if not _check_status(rkey_status, rkey_err, _tmp_files):
                logger.warn("Skipping run {}".format(run))
                continue
            logger.info("Generated randon phrase key for run {}".format(run))
            
            # Calculate md5 sum pre encryption
            if not force:
                logger.info("Calculating md5sum before encryption")
                md5_proc = sp.Popen(['md5sum', run_zip ], stdout=sp.PIPE, stderr=sp.PIPE)
                md5_status = md5_proc.wait()
                md5_out, md5_err = md5_proc.communicate()
                if not _check_status(md5_status, md5_err, _tmp_files):
                    logger.warn("Skipping run {}".format(run))
                    continue
                pre_encrypt_md5 = md5_out.split()[0]
            
            # Encrypt the zipped run file
            logger.info("Encrypting the zipped run file")
            gpg_proc = sp.Popen(['gpg', '--symmetric', '--cipher-algo', 'aes256', '--passphrase-file', run_key, '--batch', '--compress-algo', 'none', '-o', run_zip_encrypted, run_zip], 
                                stdout=sp.PIPE, stderr=sp.PIPE)
            gpg_status = gpg_proc.wait()
            gpg_out, gpg_err = gpg_proc.communicate()
            if not _check_status(gpg_status, gpg_err, _tmp_files):
                logger.warn("Skipping run {}".format(run))
                continue
            
            # Decrypt and check for md5
            if not force:
                logger.info("Calculating md5sum after encryption")
                dpg_proc = sp.Popen(['gpg', '--decrypt', '--cipher-algo', 'aes256', '--passphrase-file', run_key, '--batch', run_zip_encrypted], stdout=sp.PIPE, stderr=sp.PIPE)
                md5_proc = sp.Popen(['md5sum'], stdin=dpg_proc.stdout, stdout=sp.PIPE, stderr=sp.PIPE)
                dpg_status = dpg_proc.wait()
                md5_status = md5_proc.wait()
                md5_out, md5_err = md5_proc.communicate()
                if not _check_status(md5_status, md5_err, _tmp_files):
                    logger.warn("Skipping run {}".format(run))
                    continue
                post_encrypt_md5 = md5_out.split()[0]
                if pre_encrypt_md5 != post_encrypt_md5:
                    logger.error("md5sum did not match before ({}) and after ({}) encryption. Will remove temp files and move on".format(pre_encrypt_md5, post_encrypt_md5))
                    _clean_tmp_files(_tmp_files)
                    continue
           
            logger.info("Md5sum is macthing before and after encryption") 
            # Encrypt and move the key file
            try:
                sp.check_call(['gpg', '-e', '-r', 'Senthilkumar', '-o', run_key_encrypted, run_key])
                shutil.move(run_key_encrypted, os.path.join(keys_path, run_key_encrypted))
            except:
                logger.error("Encrption of key file failed, skipping run")
                _clean_tmp_files(_tmp_files)
                continue
            
            logger.info("Encryption of run {} is successfully done, removing zipped run file".format(run))
            _clean_tmp_files([run_zip, run_key,run_flag])
            

#############################################################
# Class helper methods, not exposed as commands/subcommands #
#############################################################

def _get_runs_from_path(path, ext=None, filter_by_ext=False):
    """Collect runs directory from given path, if an 'ext' is given
    files with those extensions are included as well. 'ext' can also
    be used to filter the result.
    
    :param str path: a path to look for runs/files
    :param str ext: an extension string for files to be included with runs
    :param boolean filter_by_ext: to use 'ext' as filter
    """
    
    runs = []
    if not os.path.exists(path) or not os.path.isdir(path):
        logger.warn("Path {} does not exist or it is not a directory".format(path))
        return runs
    
    for item in os.listdir(path):
        if filter_by_ext and ext and not item.endswith(ext):
            continue
        elif ext and item.endswith(ext):
            item = item.replace(ext, '')
        elif not os.path.isdir(item):
            continue
        if re.match(filesystem.RUN_RE, item) and item not in runs:
            runs.append(item)
    
    runs = [os.path.join(path,run) for run in runs]
    
    return runs

def _avail_disk_space(path, run, data_dirs):
    """Check the space on file system based on parent directory of the run
    
    :param str path: path to use to get the file system stat
    :param str run: run name a.k.a. FC name to get run type
    :param list data_dirs: list od dir paths to fetch ongoing runs
    """
    
    # not able to fetch runtype use the max size as precaution, size units in GB
    illumina_run_sizes = {'hiseq' : 500, 'hiseqx' : 900, 'miseq' : 20}
    required_size = illumina_run_sizes.get(_get_run_type(run), 900) * 2
    
    # check for any ongoing runs and add up the required size accrdingly
    for ddir in data_dirs:
        for item in os.listdir(ddir):
            if not re.match(filesystem.RUN_RE, item):
                continue
            if not os.path.exists(os.path.join(ddir, item, "RTAComplete.txt")):
                required_size += illumina_run_sizes.get(_get_run_type(run), 900)
    
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

def _get_run_type(run):
    """Returns run type based on the flowcell name
    
    :param str run: run name a.k.a. full FC name
    """
    
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

def _run_is_demuxed(run, couch_info):
    """Check in StatusDB 'x_flowcells' database if the given run has an entry
    which means it was demultiplexed (as TACA only creates a document upon
    successfull demultiplexing)
    
    :param str run: a run a.k.a. flowcell name
    """
    
    run_terms = run.split('_')
    run_date = run_terms[0]
    run_fc = run_terms[-1]
    run_name = "{}_{}".format(run_date, run_fc)
    
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

def _check_status(status, err_msg, files_to_remove=[]):
    """Check if a non-zero status is given and exit code
    
       :param status: returcode of Popen object
       :param file_to_remove: path of file to remove before exiting"""
    if status != 0:
        _clean_tmp_files(files_to_remove)
        logger.error("Encryption failed with the error \n{}".format(err_msg))
        return False
    return True

def _clean_tmp_files(files):
    """Remove the file is exist"""
    for fl in files:
        if os.path.exists(fl):
            os.remove(fl)
