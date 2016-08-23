"""Storage methods and utilities"""
import getpass
import logging
import os
import re
import shutil
import time

from datetime import datetime
from glob import glob
from multiprocessing import Pool

from statusdb.db import connections as statusdb
from taca.utils.config import CONFIG
from taca.utils import filesystem, misc

logger = logging.getLogger(__name__)

# This is used by many of the functions in this module
finished_run_indicator = CONFIG.get('storage', {}).get('finished_run_indicator', 'RTAComplete.txt')

def cleanup_nas(seconds):
    """Will move the finished runs in NASes to nosync directory.

    :param int seconds: Days/hours converted as second to consider a run to be old
    """
    couch_info = CONFIG.get('statusdb')
    mail_recipients = CONFIG.get('mail', {}).get('recipients')
    check_demux = CONFIG.get('storage', {}).get('check_demux', False)
    host_name = os.getenv('HOSTNAME', os.uname()[1]).split('.', 1)[0]
    for data_dir in CONFIG.get('storage').get('data_dirs'):
        logger.info('Moving old runs in {}'.format(data_dir))
        with filesystem.chdir(data_dir):
            for run in [r for r in os.listdir(data_dir) if re.match(filesystem.RUN_RE, r)]:
                rta_file = os.path.join(run, finished_run_indicator)
                if os.path.exists(rta_file):
                    if check_demux:
                        if misc.run_is_demuxed(run, couch_info):
                            logger.info('Moving run {} to nosync directory'.format(os.path.basename(run)))
                            shutil.move(run, 'nosync')
                        elif os.stat(rta_file).st_mtime < time.time() - seconds:
                            logger.warn('Run {} is older than given time, but it is not demultiplexed yet'
                                        .format(run))
                            sbt = "Run not demultiplexed - {}".format(run)
                            msg = ("Run '{}' in '{}' is older then given threshold, but seems like it is not "
                                  "yet demultiplexed".format(os.path.join(data_dir, run), host_name))
                            misc.send_mail(sbt, msg, mail_recipients)
                    else:
                        if os.stat(rta_file).st_mtime < time.time() - seconds:
                            logger.info('Moving run {} to nosync directory'.format(os.path.basename(run)))
                            shutil.move(run, 'nosync')
                        else:
                            logger.info('{} file exists but is not older than given time, skipping run {}'
                                        .format(finished_run_indicator, run))


def cleanup_processing(seconds):
    """Cleanup runs in processing server.

    :param int seconds: Days/hours converted as second to consider a run to be old
    """
    try:
        #Remove old runs from archiving dirs
        for archive_dir in CONFIG.get('storage').get('archive_dirs').values():
            logger.info('Removing old runs in {}'.format(archive_dir))
            with filesystem.chdir(archive_dir):
                for run in [r for r in os.listdir(archive_dir) if re.match(filesystem.RUN_RE, r)]:
                    rta_file = os.path.join(run, finished_run_indicator)
                    if os.path.exists(rta_file):
                        if os.stat(rta_file).st_mtime < time.time() - seconds:
                            logger.info('Removing run {} to nosync directory'.format(os.path.basename(run)))
                            shutil.rmtree(run)
                        else:
                            logger.info('{} file exists but is not older than given time, skipping run {}'.format(
                                        finished_run_indicator, run))
    except IOError:
        sbj = "Cannot archive old runs in processing server"
        msg = ("Could not find transfer.tsv file, so I cannot decide if I should "
               "archive any run or not.")
        cnt = CONFIG.get('contact', None)
        if not cnt:
            cnt = "{}@localhost".format(getpass.getuser())
        logger.error(msg)
        misc.send_mail(sbj, msg, cnt)


def cleanup_uppmax(site, seconds, dry_run=False):
    """Remove project/run that have been closed more than given time (as seconds)
    from the given 'site' on uppmax

    :param str site: site where the cleanup should be performed
    :param int seconds: Days/hours converted as second to consider a run to be old
    """
    seconds = check_default(site, seconds, CONFIG)
    if not seconds:
        return
    root_dir = CONFIG.get('cleanup').get(site).get('root')
    deleted_log = CONFIG.get('cleanup').get('deleted_log')
    assert os.path.exists(os.path.join(root_dir,deleted_log)), "Log directory {} doesn't exist in {}".format(deleted_log,root_dir)
    log_file = os.path.join(root_dir,"{fl}/{fl}.log".format(fl=deleted_log))
    list_to_delete = []

    ## get glob path patterns to search and remove from root directory
    try:
        archive_config = CONFIG['cleanup']['archive']
        ## the glob path should be relative to the run folder, like "Unaligned_*/Project_*"
        config_ppath = archive_config['proj_path']
        ## Glob path should be relative to run folder, like "Unaligned_0bp/Undetermined_indices/*/*.fastq.gz"
        config_npath = archive_config['undet_noindex']
        ## Glob path should be relative to run folder, like "Unaligned_*bp/Undetermined_indices/*/*.fastq.gz"
        config_upath = archive_config['undet_all']
    except KeyError as e:
        logger.error("Config file is missing the key {}, make sure it have all required information".format(str(e)))
        raise SystemExit

    # make a connection for project db #
    pcon = statusdb.ProjectSummaryConnection()
    assert pcon, "Could not connect to project database in StatusDB"

    if site in ["analysis", "illumina"]:
        ## work flow for cleaning up illumina/analysis ##
        projects = [ p for p in os.listdir(root_dir) if re.match(filesystem.PROJECT_RE,p) ]
        list_to_delete.extend(get_closed_projects(projects, pcon, seconds))
    elif site == "archive":
        ##work flow for cleaning archive ##
        runs = [ r for r in os.listdir(root_dir) if re.match(filesystem.RUN_RE,r) ]
        for run in runs:
            with filesystem.chdir(os.path.join(root_dir, run)):
                ## Collect all project path from demultiplexed directories in the run folder
                all_proj_path = glob(config_ppath)
                all_proj_dict = {os.path.basename(pp).replace('Project_','').replace('__', '.'): pp for pp in all_proj_path}
                closed_projects = get_closed_projects(all_proj_dict.keys(), pcon, seconds)
                ## Only proceed cleaning the data for closed projects
                for closed_proj in closed_projects:
                    closed_proj_fq = glob("{}/*/*.fastq.gz".format(all_proj_dict[closed_proj]))
                    list_to_delete.extend([os.path.join(run, pfile) for pfile in closed_proj_fq])
                ## Remove the undetermined fastq files for NoIndex case always
                undetermined_fastq_files = glob(config_npath)
                ## Remove undeterminded fastq files for all index length if all project run in the FC is closed
                if len(all_proj_dict.keys()) == len(closed_projects):
                    undetermined_fastq_files = glob(config_upath)
                list_to_delete.extend([os.path.join(run, ufile) for ufile in undetermined_fastq_files])

    ## delete and log
    for item in list_to_delete:
        if dry_run:
            logger.info('Will remove {} from {}'.format(item,root_dir))
            continue
        try:
            to_remove = os.path.join(root_dir,item)
            if os.path.isfile(to_remove):
                os.remove(to_remove)
            elif os.path.isdir(to_remove):
                shutil.rmtree(to_remove)
            logger.info('Removed {} from {}'.format(item,root_dir))
            with open(log_file,'a') as to_log:
                to_log.write("{}\t{}\n".format(to_remove,datetime.strftime(datetime.now(),'%Y-%m-%d %H:%M')))
        except OSError:
            logger.warn("Could not remove {} from {}".format(item,root_dir))
            continue


#############################################################
# Class helper methods, not exposed as commands/subcommands #
#############################################################

def get_closed_projects(projs, pj_con, seconds):
    """Takes list of project and gives project list that are closed
    more than given time(as seconds)

    :param list projs: list of projects to check
    :param obj pj_con: connection object to project database
    :param int seconds: Days/hours converted as seconds to check
    """
    closed_projs = []
    for proj in projs:
        if proj not in pj_con.name_view.keys():
            logger.warn("Project {} is not in database, so SKIPPING it.."
                        .format(proj))
            continue
        proj_db_obj = pj_con.get_entry(proj)
        try:
            proj_close_date = proj_db_obj['close_date']
        except KeyError:
            logger.warn("Project {} is either open or too old, so SKIPPING it..".format(proj))
            continue
        if misc.to_seconds(days=misc.days_old(proj_close_date,date_format='%Y-%m-%d')) > seconds:
            closed_projs.append(proj)
    return closed_projs


def check_default(site, seconds, config):
    """Check if time(as seconds) given while running command. If not take the default threshold
    from config file (which should exist). Also when 'days' given on the command line
    raise a check to make sure it was really meant to do so

    :param str site: site to be cleaned and relevent date to pick
    :param int seconds: Days/hours converted as seconds to check
    :param dict config: config file parsed and saved as dictionary
    """
    try:
        default_days = config['cleanup'][site]['days']
        default_seconds = misc.to_seconds(days=default_days)
    except KeyError:
        raise
    if not seconds:
        return default_seconds
    elif seconds >= default_seconds:
        return seconds
    else:
        if misc.query_yes_no("Seems like given time is less than the "
                             " default({}) days, are you sure to proceed ?"
                             .format(default_days), default="no"):
            return seconds
        else:
            return None


############################################################
#######             DEPRECATED METHODS               #######
############################################################
## these methods are not in use anymore but kept here as it
## might be used again in future, when world is about to end

def archive_to_swestore(seconds, run=None, max_runs=None, force=False, compress_only=False):
    """Send runs (as archives) in NAS nosync to swestore for backup

    :param int seconds: Days/hours converted as seconds to check
    :param str run: specific run to send swestore
    :param int max_runs: number of runs to be processed simultaneously
    :param bool force: Force the archiving even if the run is not complete
    :param bool compress_only: Compress the run without sending it to swestore
    """
    # If the run is specified in the command line, check that exists and archive
    if run:
        run = os.path.basename(run)
        base_dir = os.path.dirname(run)
        if re.match(filesystem.RUN_RE, run):
            # If the parameter is not an absolute path, find the run in the archive_dirs
            if not base_dir:
                for archive_dir in CONFIG.get('storage').get('archive_dirs'):
                    if os.path.exists(os.path.join(archive_dir, run)):
                        base_dir = archive_dir
            if not os.path.exists(os.path.join(base_dir, run)):
                logger.error(("Run {} not found. Please make sure to specify "
                              "the absolute path or relative path being in "
                              "the correct directory.".format(run)))
            else:
                with filesystem.chdir(base_dir):
                    _archive_run((run, seconds, force, compress_only))
        else:
            logger.error("The name {} doesn't look like an Illumina run"
                         .format(os.path.basename(run)))
    # Otherwise find all runs in every data dir on the nosync partition
    else:
        logger.info("Archiving old runs to SWESTORE")
        for to_send_dir in CONFIG.get('storage').get('archive_dirs'):
            logger.info('Checking {} directory'.format(to_send_dir))
            with filesystem.chdir(to_send_dir):
                to_be_archived = [r for r in os.listdir(to_send_dir)
                                  if re.match(filesystem.RUN_RE, r)
                                  and not os.path.exists("{}.archiving".format(r.split('.')[0]))]
                if to_be_archived:
                    pool = Pool(processes=len(to_be_archived) if not max_runs else max_runs)
                    pool.map_async(_archive_run, ((run, seconds, force, compress_only) for run in to_be_archived))
                    pool.close()
                    pool.join()
                else:
                    logger.info('No old runs to be archived')


def _archive_run((run, seconds, force, compress_only)):
    """ Archive a specific run to swestore

    :param str run: Run directory
    :param int seconds: Days/hours converted as seconds to check
    :param bool force: Force the archiving even if the run is not complete
    :param bool compress_only: Only compress the run without sending it to swestore
    """

    def _send_to_swestore(f, dest, remove=True):
        """ Send file to swestore checking adler32 on destination and eventually
        removing the file from disk

        :param str f: File to remove
        :param str dest: Destination directory in Swestore
        :param bool remove: If True, remove original file from source
        """
        if not filesystem.is_in_swestore(f):
            logger.info("Sending {} to swestore".format(f))
            misc.call_external_command('iput -R swestoreArchCacheResc -P {file} {dest}'.format(file=f, dest=dest),
                    with_log_files=True, prefix=f.replace('.tar.bz2',''), log_dir="swestore_logs")
            logger.info('Run {} sent to swestore.'.format(f))
            if remove:
                logger.info('Removing run'.format(f))
                os.remove(f)
        else:
            logger.warn('Run {} is already in Swestore, not sending it again nor removing from the disk'.format(f))

    # Create state file to say that the run is being archived
    open("{}.archiving".format(run.split('.')[0]), 'w').close()
    if run.endswith('bz2'):
        if os.stat(run).st_mtime < time.time() - seconds:
            _send_to_swestore(run, CONFIG.get('storage').get('irods').get('irodsHome'))
        else:
            logger.info("Run {} is not older than given time yet. Not archiving".format(run))
    else:
        rta_file = os.path.join(run, finished_run_indicator)
        if not os.path.exists(rta_file) and not force:
            logger.warn(("Run {} doesn't seem to be completed and --force option was "
                      "not enabled, not archiving the run".format(run)))
        if force or (os.path.exists(rta_file) and os.stat(rta_file).st_mtime < time.time() - seconds):
            logger.info("Compressing run {}".format(run))
            # Compress with pbzip2
            misc.call_external_command('tar --use-compress-program=pbzip2 -cf {run}.tar.bz2 {run}'.format(run=run))
            logger.info('Run {} successfully compressed! Removing from disk...'.format(run))
            shutil.rmtree(run)
            if not compress_only:
                _send_to_swestore('{}.tar.bz2'.format(run), CONFIG.get('storage').get('irods').get('irodsHome'))
        else:
            logger.info("Run {} is not completed or is not older than given time yet. Not archiving".format(run))
    os.remove("{}.archiving".format(run.split('.')[0]))


def cleanup_swestore(seconds, dry_run=False):
    """Remove archived runs from swestore

    :param int seconds: Days/hours converted as seconds to check
    """
    seconds = check_default(site, seconds, CONFIG)
    if not seconds:
        return
    runs = filesystem.list_runs_in_swestore(path=CONFIG.get('cleanup').get('swestore').get('root'))
    for run in runs:
        date = run.split('_')[0]
        if misc.to_seconds(misc.days_old(date)) > seconds:
            if dry_run:
                logger.info('Will remove file {} from swestore'.format(run))
                continue
            misc.call_external_command('irm -f {}'.format(run))
            logger.info('Removed file {} from swestore'.format(run))
