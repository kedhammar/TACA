"""Storage methods and utilities"""
import getpass
import logging
import os
import re
import shutil
import time

from collections import defaultdict
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


def cleanup_milou(site, seconds, dry_run=False):
    """Remove project/run that have been closed more than given time (as seconds)
    from the given 'site' on uppmax

    :param str site: site where the cleanup should be performed
    :param int seconds: Days/hours converted as second to consider a run to be old
    :param bool dry_run: Will summarize what is going to be done without really doing it
    """
    seconds = check_default(site, seconds, CONFIG)
    if not seconds:
        return
    root_dir = CONFIG.get('cleanup').get('milou').get(site).get('root')
    deleted_log = CONFIG.get('cleanup').get('milou').get('deleted_log')
    assert os.path.exists(os.path.join(root_dir,deleted_log)), "Log directory {} doesn't exist in {}".format(deleted_log,root_dir)
    log_file = os.path.join(root_dir,"{fl}/{fl}.log".format(fl=deleted_log))
    list_to_delete = []

    ## get glob path patterns to search and remove from root directory
    try:
        archive_config = CONFIG['cleanup']['milou']['archive']
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


def cleanup_irma(conf_file, days_fastq, days_analysis, only_fastq, only_analysis, dry_run=False):
    """Remove fastq/analysis data for projects that have been closed more than given 
    days (as days_fastq/days_analysis) from the given 'irma' cluster

    :param int days_fastq: Days to consider to remove fastq files for project
    :param int days_analysis: Days to consider to remove analysis data for project
    :param bool only_fastq: Remove only fastq files for closed projects
    :param bool only_analysis: Remove only analysis data for closed projects
    :param bool dry_run: Will summarize what is going to be done without really doing it
    
    Example for mat for config file
    cleanup:
        irma:
            flowcell:
                ##this path is nothing but incoming directory, can given multiple paths
                root: 
                    - path/to/flowcells_dir
                relative_project_source: Demultiplexing
    
            ##this is path where projects are organized
            data_dir: path/to/data_dir
            analysis:
                ##directory where analysis are perfoemed for projects
                root: path/to/analysis_dir
                #should be exactly same as the qc folder name and files wished to be removed
                files_to_remove:
                    piper_ngi: 
                        - "*.bam"
    """
    try:
        config = CONFIG['cleanup']['irma']
        flowcell_dir_root = config['flowcell']['root']
        flowcell_project_source = config['flowcell']['relative_project_source']
        data_dir = config['data_dir']
        analysis_dir = config['analysis']['root']
        analysis_data_to_remove = config['analysis']['files_to_remove']
    except KeyError as e:
        logger.error("Config file is missing the key {}, make sure it have all required information".format(str(e)))
        raise SystemExit
    
    # make a connection for project db #
    pcon = statusdb.ProjectSummaryConnection(conf=conf_file)
    assert pcon, "Could not connect to project database in StatusDB"

    #compile list for project to delete
    project_clean_list, project_processed_list = ({}, [])
    logger.info("Building initial project list for removing data..")
    if only_fastq:
        logger.info("Option 'only_fastq' is given, so will not look for analysis data")
    elif only_analysis:
        logger.info("Option 'only_analysis' is given, so will not look for fastq data")
     
    if only_analysis:
        for pid in [d for d in os.listdir(analysis_dir) if re.match(r'^P\d+$', d) and \
                    not os.path.exists(os.path.join(analysis_dir, d, "cleaned"))]:
            proj_abs_path = os.path.join(analysis_dir, pid)
            proj_info = get_closed_proj_info(pcon.get_entry(pid, use_id_view=True))
            if proj_info and proj_info['closed_days'] >= days_analysis:
                analysis_data, analysis_size = collect_analysis_data_irma(pid, analysis_dir, analysis_data_to_remove)
                proj_info['analysis_to_remove'] = analysis_data
                proj_info['analysis_size'] = analysis_size
                proj_info['fastq_to_remove'] = "not_selected"
                proj_info['fastq_size'] = 0
                project_clean_list[proj_info['name']] = proj_info
    else:
        for flowcell_dir in flowcell_dir_root:
            for fc in [d for d in os.listdir(flowcell_dir) if re.match(filesystem.RUN_RE,d)]:
                fc_abs_path = os.path.join(flowcell_dir, fc)
                with filesystem.chdir(fc_abs_path):
                    projects_in_fc = [d for d in os.listdir(flowcell_project_source) \
                                      if re.match(r'^[A-Z]+[_\.]+[A-Za-z]+_\d\d_\d\d$',d) and \
                                      not os.path.exists(os.path.join(flowcell_project_source, d, "cleaned"))]
                    for _proj in projects_in_fc:
                        proj = re.sub(r'_+', '.', _proj, 1)
                        # if a project is already processed no need of fetching it again from status db
                        if proj in project_processed_list:
                            # if the project is closed more than threshold days collect the fastq files from FC
                            # no need of looking for analysis data as they would have been collected in the first time
                            if proj in project_clean_list and project_clean_list[proj]['closed_days'] >= days_fastq:
                                fc_fq_files, fq_size = collect_fastq_data_irma(fc_abs_path, os.path.join(flowcell_project_source, _proj))
                                project_clean_list[proj]['fastq_to_remove']['flowcells'][fc] = fc_fq_files['flowcells'][fc]
                                project_clean_list[proj]['fastq_size'] += fq_size
                            continue
                        project_processed_list.append(proj)
                        #by default assume all projects are not old enough for delete
                        fastq_data, analysis_data = ("young", "young")
                        fastq_size, analysis_size = (0, 0)
                        proj_info = get_closed_proj_info(pcon.get_entry(proj))
                        if proj_info:
                            # if project not old enough for fastq files and only fastq files selected move on to next project
                            if proj_info['closed_days'] >= days_fastq:
                                fastq_data, fastq_size = collect_fastq_data_irma(fc_abs_path, os.path.join(flowcell_project_source, _proj),
                                                                                 data_dir, proj_info['pid'])
                            if not only_fastq:
                                # if project is old enough for fastq files and not 'only_fastq' try collect analysis files 
                                if proj_info['closed_days'] >= days_analysis:
                                    analysis_data, analysis_size = collect_analysis_data_irma(proj_info['pid'], analysis_dir, analysis_data_to_remove)
                                # if both fastq and analysis files are not old enough move on
                                if (analysis_data == fastq_data) or ((not analysis_data or analysis_data == "cleaned") and fastq_data == "young"):
                                    continue
                            elif fastq_data == "young":
                                continue
                            else:
                                analysis_data = "not_selected"
                            proj_info['fastq_to_remove'] = fastq_data
                            proj_info['fastq_size'] = fastq_size
                            proj_info['analysis_to_remove'] = analysis_data
                            proj_info['analysis_size'] = analysis_size
                            project_clean_list[proj] = proj_info
    
    if not project_clean_list:
        logger.info("There are no projects to clean")
        return
                    
    get_files_size_text(project_clean_list)
    logger.info("Initial list is built with {} projects {}".format(len(project_clean_list), get_files_size_text(project_clean_list)))
    if  misc.query_yes_no("Interactively filter projects for cleanup ?", default="yes"):
        filtered_project = []
        #go through complied project list and remove files
        for proj, info in project_clean_list.iteritems():
            if not misc.query_yes_no("{}Delete files for this project".format(get_proj_meta_info(info, days_fastq)), default="no"):
                logger.info("Will not remove files for project {}".format(proj))
                filtered_project.append(proj)
        # remove projects that were decided not to delete
        map(project_clean_list.pop, filtered_project)
        logger.info("Removed {} projects from initial list".format(len(filtered_project)))
        logger.info("Final list is created with {} projects {}".format(len(project_clean_list), get_files_size_text(project_clean_list)))
        if not misc.query_yes_no("Proceed with cleanup ?", default="no"):
            logger.info("Aborting cleanup")
            return
    logger.info("Will start cleaning up project now")
    
    for proj, info in project_clean_list.iteritems():
        fastq_info = info.get('fastq_to_remove')
        if fastq_info and isinstance(fastq_info, dict):
            logger.info("Cleaning fastq files for project {}".format(proj))
            fastq_fc = fastq_info.get('flowcells', {})
            for fc, fc_info in fastq_fc.iteritems():
                proj_fc_root = fc_info['proj_root']
                logger.info("Removing fastq files from {}".format(proj_fc_root))
                if not dry_run:
                    _remove_files(fc_info['fq_files'])
                    logger.info("Removed fastq files from FC {} for project {}, marking it as cleaned".format(fc, proj))
                    _touch_cleaned(proj_fc_root)
            fastq_proj = fastq_info.get('proj_data')
            if fastq_proj and isinstance(fastq_proj, dict):
                proj_data_root = fastq_proj['proj_data_root']
                logger.info("Removing fastq_files from {}".format(proj_data_root))
                if not dry_run:
                    _remove_files(fastq_proj['fastq_files'])
                    logger.info("Removed fastq files from data directory for project {}, marking it as cleaned".format(proj))
                    _touch_cleaned(proj_data_root)
            
        analysis_info = info.get('analysis_to_remove')
        if analysis_info and isinstance(analysis_info, dict):
            proj_analysis_root = analysis_info['proj_analysis_root']
            logger.info("cleaning analysis data for project {}".format(proj))
            for qc, files in analysis_info['analysis_files'].iteritems():
                removed_qc = []
                logger.info("Removing files of '{}' from {}".format(qc, proj_analysis_root))
                if not dry_run:
                    _remove_files(files)
                    removed_qc.append(qc)
            map(analysis_info['analysis_files'].pop, removed_qc)
            if len(analysis_info['analysis_files']) == 0:
                logger.info("Removed analysis data for project {}, marking it cleaned".format(proj))
                _touch_cleaned(proj_analysis_root)


#############################################################
# Class helper methods, not exposed as commands/subcommands #
#############################################################

def get_closed_proj_info(pdoc):
    """check and return a dict if project is closed"""
    pdict = None
    if "close_date" in pdoc:
        closed_date = pdoc['close_date']
        closed_days = misc.days_old(closed_date, "%Y-%m-%d")
        if closed_days and isinstance(closed_days, int):
            pdict = {'name' : pdoc.get('project_name'),
                     'pid' : pdoc.get('project_id'),
                     'closed_date' : closed_date,
                     'closed_days' : closed_days,
                     'bioinfo_responsible' : pdoc.get('project_summary',{}).get('bioinfo_responsible','')}
        else:
            logger.warn("Problem calculating closed days for project {} with close data {}. Skipping it".format(
                        pdoc.get('project_name'), closed_date))
                     
    return pdict

def collect_analysis_data_irma(pid, analysis_root, files_ext_to_remove={}):
    """Collect the analysis files that have to be removed from IRMA
    return a tuple with files and total size of collected files"""
    size = 0
    proj_abs_path = os.path.join(analysis_root, pid)
    if not os.path.exists(proj_abs_path):
        file_list = None
    elif os.path.exists(os.path.join(proj_abs_path, "cleaned")):
        file_list = "cleaned"
    else:
        file_list = {'proj_analysis_root':proj_abs_path,
                     'analysis_files': defaultdict(list)}
        for qc_type,ext in files_ext_to_remove.items():
            qc_path = os.path.join(proj_abs_path, qc_type)
            if os.path.exists(qc_path):
                file_list['analysis_files'][qc_type].extend(collect_files_by_ext(qc_path, ext))
    try:
        size += sum([sum(map(os.path.getsize, fls)) for fls in file_list['analysis_files'].values()])
    except:
        pass
    return (file_list, size)

def collect_fastq_data_irma(fc_root, fc_proj_src, proj_root=None, pid=None):
    """Collect the fastq files that have to be removed from IRMA
    return a tuple with files and total size of collected files"""
    size = 0
    file_list = {'flowcells': defaultdict(dict)}
    fc_proj_path = os.path.join(fc_root, fc_proj_src)
    fc_id = os.path.basename(fc_root)
    file_list['flowcells'][fc_id] = {'proj_root': fc_proj_path,
                                     'fq_files': collect_files_by_ext(fc_proj_path, "*.fastq.gz")}
    if proj_root and pid:
        proj_abs_path = os.path.join(proj_root, pid)
        if not os.path.exists(proj_abs_path):
            file_list['proj_data'] = None
        elif os.path.exists(os.path.join(proj_abs_path, "cleaned")):
            file_list['proj_data'] = "cleaned"
        else:
            file_list['proj_data'] = {'proj_data_root': proj_abs_path,
                                      'fastq_files' : collect_files_by_ext(proj_abs_path, "*.fastq.gz")}
    size += sum(map(os.path.getsize, file_list['flowcells'][fc_id]['fq_files']))
    try:
        size += sum(map(os.path.getsize, file_list['proj_data']['fastq_files']))
    except:
        pass
    return (file_list, size)

def collect_files_by_ext(path, ext=[]):
    """Collect files with given extension from given path"""
    if isinstance(ext, str):
        ext = [ext]
    collected_files = []
    for root, dirs, files in os.walk(path):
        for e in ext:
            collected_files.extend(glob(os.path.join(root,e)))
        for d in dirs:
            collected_files.extend(collect_files_by_ext(d, ext))
    return collected_files

def get_proj_meta_info(info, days_fastq):
    """From given info collect meta info for a project"""    
    template = ("\nProject overview - {proj_name}\n"
                "Project id: {proj_id}\n"
                "Bioinfo Responsible: {bio_res}\n"
                "Closed for (days): {closed_days}\n"
                "Closed from (date): {closed_date}\n")
    
    meta_info = {'proj_name' : info.get('name'),
                 'proj_id' : info.get('pid'),
                 'bio_res' : info.get('bioinfo_responsible'),
                 'closed_days' : info.get('closed_days'),
                 'closed_date' : info.get('closed_date')}
    
    template = template.format(**meta_info)
    
    # set analysis info based upon what we have
    analysis_info = info.get('analysis_to_remove')
    if not analysis_info:
        template += "Project analysis: No analysis directory\n"
    elif isinstance(analysis_info, str) and analysis_info == "cleaned":
        template += "Project analysis: Analysis directory already cleaned\n"
    elif isinstance(analysis_info, dict):
        f_stat = []
        for qc_type, files in analysis_info['analysis_files'].iteritems():
            f_stat.append("{} ({} files)".format(qc_type, len(files)))
        template += "Project analyzed: {}\n".format(", ".join(f_stat))
    
    # set fastq info based upon what we have
    fq_info = info.get('fastq_to_remove')
    if isinstance(fq_info, str) and fq_info == "young":
        template += "Project been closed less than {} days, so will not remove any 'fastq' files\n".format(days_fastq)
    elif isinstance(fq_info, dict):
        proj_fq_info = fq_info.get('proj_data')
        if not proj_fq_info:
            template += "Project organized: No organized directory for project\n"
        elif isinstance(proj_fq_info, str) and proj_fq_info == "cleaned":
            template += "Project organized: Project directory is already cleaned\n"
        elif isinstance(proj_fq_info, dict):
            template += "Project organized: Project is organized with {} fastq files\n".format(len(proj_fq_info['fastq_files']))
        fc_fq_info = fq_info.get('flowcells', {})
        fc_num = len(fc_fq_info.keys())
        fc_files = sum(map(len, [fc_info.get('fq_files', [])for fc_info in fc_fq_info.values()]))
        template += "Flowcells: There are {} FC with total {} fastq files\n".format(fc_num, fc_files)

    return template

def get_files_size_text(plist):
    """Get project list dict and give back string with overll sizes"""
    def _def_get_size_unit(s):
        kb = 1000
        mb = kb * 1000
        gb = mb * 1000
        tb = gb * 1000
        if s > tb:
            s = "~{}tb".format(s/tb)
        elif s > gb:
            s = "~{}gb".format(s/gb)
        elif s > mb:
            s = "~{}mb".format(s/mb)
        elif s > kb:
            s = "~{}kb".format(s/kb)
        elif s > 0:
            s = "~{}b".format(s/b)
        return s
    fsize = _def_get_size_unit(sum([i.get('fastq_size',0) for i in plist.values()]))
    asize = _def_get_size_unit(sum([i.get('analysis_size',0) for i in plist.values()]))
    return "({f}{s}{a}) ".format(f = "~{} fastq data".format(fsize) if fsize else "",
                                 a = "~{} analysis data".format(asize) if asize else "",
                                 s = " and " if fsize and asize else "")

def _remove_files(files):
    """Remove files from given list"""
    for fl in files:
        os.remove(fl)

def _touch_cleaned(path):
    """Touch a 'cleaned' file in a given path"""
    open(os.path.join(path, "cleaned"), 'w').close()

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
        default_days = config['cleanup']['milou'][site]['days']
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
