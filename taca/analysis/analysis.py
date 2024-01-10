"""Analysis methods for TACA."""
import glob
import logging
import os
import sys
import subprocess

from shutil import copyfile
from shutil import copytree
from taca.illumina.Standard_Runs import Standard_Run
from taca.illumina.MiSeq_Runs import MiSeq_Run
from taca.illumina.NextSeq_Runs import NextSeq_Run
from taca.illumina.NovaSeq_Runs import NovaSeq_Run
from taca.illumina.NovaSeqXPlus_Runs import NovaSeqXPlus_Run
from taca.utils.config import CONFIG
from taca.utils.transfer import RsyncAgent
from taca.utils import statusdb

from flowcell_parser.classes import RunParametersParser
from io import open

logger = logging.getLogger(__name__)


def get_runObj(run, software):
    """Tries to read runParameters.xml to parse the type of sequencer
        and then return the respective Run object (MiSeq, HiSeq..)

    :param run: run name identifier
    :type run: string
    :returns: returns the sequencer type object,
    None if the sequencer type is unknown of there was an error
    """

    if os.path.exists(os.path.join(run, 'runParameters.xml')):
        run_parameters_file = 'runParameters.xml'
    elif os.path.exists(os.path.join(run, 'RunParameters.xml')):
        run_parameters_file = 'RunParameters.xml'
    else:
        logger.error('Cannot find RunParameters.xml or runParameters.xml in the run folder for run {}'.format(run))
        return

    run_parameters_path = os.path.join(run, run_parameters_file)
    try:
        run_parameters = RunParametersParser(run_parameters_path)
    except OSError:
        logger.warn('Problems parsing the runParameters.xml file at {}. '
                    'This is quite unexpected. please archive the run {} manually'.format(run_parameters_path, run))
    else:
        # Do a case by case test because there are so many version of RunParameters that there is no real other way
        runtype = run_parameters.data['RunParameters'].get('InstrumentType',
                                               run_parameters.data['RunParameters'].get('ApplicationName',
                                                                            run_parameters.data['RunParameters'].get('Application',
                                                                                                         '')))
        if 'Setup' in run_parameters.data['RunParameters']:
            # This is the HiSeq2500, MiSeq, and HiSeqX case
            try:
                # Works for recent control software
                runtype = run_parameters.data['RunParameters']['Setup']['Flowcell']
            except KeyError:
                # Use this as second resource but print a warning in the logs
                logger.warn('Parsing runParameters to fecth instrument type, '
                            'not found Flowcell information in it. Using ApplicationName')
                # Here makes sense to use get with default value '' ->
                # so that it doesn't raise an exception in the next lines
                # (in case ApplicationName is not found, get returns None)
                runtype = run_parameters.data['RunParameters']['Setup'].get('ApplicationName', '')

        if 'MiSeq' in runtype:
            return MiSeq_Run(run, software, CONFIG['analysis']['MiSeq'])
        elif 'NextSeq' in runtype:
            return NextSeq_Run(run, software, CONFIG['analysis']['NextSeq'])
        elif 'NovaSeqXPlus' in runtype:
            return NovaSeqXPlus_Run(run, software, CONFIG['analysis']['NovaSeqXPlus'])
        elif 'NovaSeq' in runtype:
            return NovaSeq_Run(run, software, CONFIG['analysis']['NovaSeq'])
        else:
            logger.warn('Unrecognized run type {}, cannot archive the run {}. '
                        'Someone as likely bought a new sequencer without telling '
                        'it to the bioinfo team'.format(runtype, run))
    return None

def upload_to_statusdb(run_dir, software):
    """Function to upload run_dir informations to statusDB directly from click interface.

    :param run_dir: run name identifier
    :type run: string
    :rtype: None
    """
    runObj = get_runObj(run_dir, software)
    if runObj:
        # runObj can be None
        # Make the actual upload
        _upload_to_statusdb(runObj)

def _upload_to_statusdb(run):
    """Triggers the upload to statusdb using the dependency flowcell_parser.

    :param Run run: the object run
    """
    couch_conf = CONFIG['statusdb']
    couch_connection = statusdb.StatusdbSession(couch_conf).connection
    db = couch_connection[couch_conf['xten_db']]
    parser = run.runParserObj
    # Check if I have NoIndex lanes
    for element in parser.obj['samplesheet_csv']:
        if 'NoIndex' in element['index'] or not element['index']: # NoIndex in the case of HiSeq, empty in the case of HiSeqX
            lane = element['Lane'] # This is a lane with NoIndex
            # In this case PF Cluster is the number of undetermined reads
            try:
                PFclusters = parser.obj['Undetermined'][lane]['unknown']
            except KeyError:
                logger.error('While taking extra care of lane {} of NoIndex type ' \
                             'I found out that not all values were available'.format(lane))
                continue
            # In Lanes_stats fix the lane yield
            parser.obj['illumina']['Demultiplex_Stats']['Lanes_stats'][int(lane) - 1]['PF Clusters'] = str(PFclusters)
            # Now fix Barcode lane stats
            updated = 0 # Check that only one update is made
            for sample in parser.obj['illumina']['Demultiplex_Stats']['Barcode_lane_statistics']:
                if lane in sample['Lane']:
                    updated += 1
                    sample['PF Clusters'] = str(PFclusters)
            if updated != 1:
                logger.error('While taking extra care of lane {} of NoIndex type '
                             'I updated more than once the barcode_lane. '
                             'This is too much to continue so I will fail.'.format(lane))
                os.sys.exit()
            # If I am here it means I changed the HTML representation to something
            # else to accomodate the wired things we do
            # someone told me that in such cases it is better to put a place holder for this
            parser.obj['illumina']['Demultiplex_Stats']['NotOriginal'] = 'True'
    # Update info about bcl2fastq tool
    if not parser.obj.get('DemultiplexConfig'):
        parser.obj['DemultiplexConfig'] = {'Setup': {'Software': run.CONFIG.get('bcl2fastq', {})}}
    statusdb.update_doc(db, parser.obj, over_write_db_entry=True)

def transfer_run(run_dir):
    """Interface for click to force a transfer a run to uppmax.

    :param: string run_dir: the run to tranfer
    """
    runObj = get_runObj(run_dir, software)
    mail_recipients = CONFIG.get('mail', {}).get('recipients')
    if runObj is None:
        mail_recipients = CONFIG.get('mail', {}).get('recipients')
        logger.error('Trying to force a transfer of run {} but the sequencer was not recognized.'.format(run_dir))
    else:
        runObj.transfer_run(os.path.join('nosync', CONFIG['analysis']['status_dir'], 'transfer.tsv'), mail_recipients)

def transfer_runfolder(run_dir, pid, exclude_lane):
    """Transfer the entire run folder for a specified project and run to uppmax.

    :param: string run_dir: the run to transfer
    :param: string pid: the projects to include in the SampleSheet separated by comma
    :param: string exclude_lane: lanes to exclude separated by comma

    """
    # Validate whether run_dir exists or is valid
    run_dir = os.path.abspath(run_dir)
    if not os.path.exists(run_dir) or not os.path.isdir(run_dir):
        logger.error('Unable to locate the specified run directory for transfer.')
        sys.exit()

    original_sample_sheet = os.path.join(run_dir, 'SampleSheet.csv')
    pid_list = list(set([x.strip() for x in pid.split(',')]))
    new_sample_sheet = os.path.join(run_dir, '_'.join(pid_list) + '_SampleSheet.txt')

    # Write new sample sheet including only rows for the specified project
    try:
        with open(new_sample_sheet, 'w') as nss:
            nss.write(extract_project_samplesheet(original_sample_sheet, pid_list))
    except IOError as e:
        logger.error('An error occured while parsing the samplesheet. '
        'Please check the sample sheet and try again.')
        raise e

    # Create a tar archive of the runfolder
    dir_name = os.path.basename(run_dir)
    archive = run_dir + '.tar.gz'
    run_dir_path = os.path.dirname(run_dir)

    # Prepare the options for excluding lanes
    if exclude_lane != '':
        dir_for_excluding_lane = []
        lane_to_exclude = exclude_lane.split(',')
        for lane in lane_to_exclude:
            if os.path.isdir('{}/{}/Thumbnail_Images/L00{}'.format(run_dir_path, dir_name, lane)):
                dir_for_excluding_lane.extend(['--exclude', 'Thumbnail_Images/L00{}'.format(lane)])
            if os.path.isdir('{}/{}/Images/Focus/L00{}'.format(run_dir_path, dir_name, lane)):
                dir_for_excluding_lane.extend(['--exclude', 'Images/Focus/L00{}'.format(lane)])
            if os.path.isdir('{}/{}/Data/Intensities/L00{}'.format(run_dir_path, dir_name, lane)):
                dir_for_excluding_lane.extend(['--exclude', 'Data/Intensities/L00{}'.format(lane)])
            if os.path.isdir('{}/{}/Data/Intensities/BaseCalls/L00{}'.format(run_dir_path, dir_name, lane)):
                dir_for_excluding_lane.extend(['--exclude', 'Data/Intensities/BaseCalls/L00{}'.format(lane)])

    try:
        exclude_options_for_tar = ['--exclude', 'Demultiplexing*',
                         '--exclude', 'demux_*',
                         '--exclude', 'rsync*',
                         '--exclude', '*.csv']
        if exclude_lane != '':
            exclude_options_for_tar += dir_for_excluding_lane

        subprocess.call(['tar'] + exclude_options_for_tar + ['-cvzf', archive, '-C', run_dir_path, dir_name])
    except subprocess.CalledProcessError as e:
        logger.error('Error creating tar archive')
        raise e

    # Generate the md5sum under the same folder as run_dir
    md5file = archive + '.md5'
    try:
        f = open(md5file, 'w')
        os.chdir(run_dir_path)
        subprocess.call(['md5sum', os.path.basename(archive)], stdout=f)
        f.close()
    except subprocess.CalledProcessError as e:
        logger.error('Error creating md5 file')
        raise e

    # Rsync the files to the analysis cluster
    destination = CONFIG['analysis']['deliver_runfolder'].get('destination')
    rsync_opts = {'-LtDrv': None,
                  '--chmod': 'g+rw'}
    connection_details = CONFIG['analysis']['deliver_runfolder'].get('analysis_server')
    archive_transfer = RsyncAgent(archive,
                                  dest_path=destination,
                                  remote_host=connection_details['host'],
                                  remote_user=connection_details['user'],
                                  validate=False,
                                  opts=rsync_opts)
    md5_transfer = RsyncAgent(md5file,
                              dest_path=destination,
                              remote_host=connection_details['host'],
                              remote_user=connection_details['user'],
                              validate=False,
                              opts=rsync_opts)

    archive_transfer.transfer()
    md5_transfer.transfer()

    # clean up the generated files
    try:
        os.remove(new_sample_sheet)
        os.remove(archive)
        os.remove(md5file)
    except IOError as e:
        logger.error('Was not able to delete all temporary files')
        raise e
    return

def extract_project_samplesheet(sample_sheet, pid_list):
    header_line = ''
    project_entries = ''
    with open(sample_sheet) as f:
        for line in f:
            if line.split(',')[0] in ('Lane', 'FCID'):  # include the header
                header_line += line
            elif any(pid in line for pid in pid_list):
                project_entries += line  # include only lines related to the specified project
    new_samplesheet_content = header_line + project_entries
    return new_samplesheet_content

def run_preprocessing(run, software):
    """Run demultiplexing in all data directories.

    :param str run: Process a particular run instead of looking for runs
    """
    def _process(run):
        """Process a run/flowcell and transfer to analysis server.

        :param taca.illumina.Run run: Run to be processed and transferred
        """
        logger.info('Checking run {}'.format(run.id))
        transfer_file = os.path.join(CONFIG['analysis']['status_dir'], 'transfer.tsv')
        if run.is_transferred(transfer_file):  # Transfer is ongoing or finished. Do nothing. Sometimes caused by runs that are copied back from NAS after a reboot
            logger.info('Run {} already transferred to analysis server, skipping it'.format(run.id))
            return

        if run.get_run_status() == 'SEQUENCING':
            logger.info('Run {} is not finished yet'.format(run.id))
            if 'statusdb' in CONFIG:
                _upload_to_statusdb(run)
        elif run.get_run_status() == 'TO_START':
            if run.get_run_type() == 'NON-NGI-RUN':
                # For now MiSeq specific case. Process only NGI-run, skip all the others (PhD student runs)
                logger.warn('Run {} marked as {}, '
                            'TACA will skip this and move the run to '
                            'no-sync directory'.format(run.id, run.get_run_type()))
                if 'storage' in CONFIG:
                    run.archive_run(CONFIG['storage']['archive_dirs'][run.sequencer_type])
                return
            logger.info(('Starting BCL to FASTQ conversion and demultiplexing for run {}'.format(run.id)))
            if 'statusdb' in CONFIG:
                _upload_to_statusdb(run)
            run.demultiplex_run()
        elif run.get_run_status() == 'IN_PROGRESS':
            logger.info(('BCL conversion and demultiplexing process in '
                         'progress for run {}, skipping it'.format(run.id)))
            # Upload to statusDB if applies
            if 'statusdb' in CONFIG:
                _upload_to_statusdb(run)
            # This function checks if demux is done
            run.check_run_status()

        # Previous elif might change the status to COMPLETED, therefore to avoid skipping
        # a cycle take the last if out of the elif
        if run.get_run_status() == 'COMPLETED':
            run.check_run_status()
            logger.info(('Preprocessing of run {} is finished, transferring it'.format(run.id)))
            # Upload to statusDB if applies
            if 'statusdb' in CONFIG:
                _upload_to_statusdb(run)
                demux_summary_message = []
                for demux_id, demux_log in run.demux_summary.items():
                    if demux_log['errors'] or demux_log['warnings']:
                        demux_summary_message.append("Sub-Demultiplexing in Demultiplexing_{} completed with {} errors and {} warnings:".format(demux_id, demux_log['errors'], demux_log['warnings']))
                        demux_summary_message.append("\n".join(demux_log['error_and_warning_messages'][:5]))
                        if len(demux_log['error_and_warning_messages'])>5:
                            demux_summary_message.append("...... Only the first 5 errors or warnings are displayed for Demultiplexing_{}.".format(demux_id))
                # Notify with a mail run completion and stats uploaded
                if demux_summary_message:
                    sbt = ("{} Demultiplexing Completed with ERRORs or WARNINGS!".format(run.id))
                    msg = """The run {run} has been demultiplexed with errors or warnings!

                    {errors_warnings}

                    The Run will be transferred to the analysis cluster for further analysis.

                    The run is available at : https://genomics-status.scilifelab.se/flowcells/{run}

                    """.format(errors_warnings='\n'.join(demux_summary_message), run=run.id)
                else:
                    sbt = ("{} Demultiplexing Completed!".format(run.id))
                    msg = """The run {run} has been demultiplexed without any error or warning.

                    The Run will be transferred to the analysis cluster for further analysis.

                    The run is available at : https://genomics-status.scilifelab.se/flowcells/{run}

                    """.format(run=run.id)
                run.send_mail(sbt, msg, rcp=CONFIG['mail']['recipients'])

            # Copy demultiplex stats file, InterOp meta data and run xml files to shared file system for LIMS purpose
            if 'mfs_path' in CONFIG['analysis']:
                try:
                    mfs_dest = os.path.join(CONFIG['analysis']['mfs_path'][run.sequencer_type.lower()],run.id)
                    logger.info('Copying demultiplex stats, InterOp metadata and XML files for run {} to {}'.format(run.id, mfs_dest))
                    if not os.path.exists(mfs_dest):
                        os.mkdir(mfs_dest)
                    demulti_stat_src = os.path.join(run.run_dir, run.demux_dir, 'Reports',
                                                    'html', run.flowcell_id, 'all', 'all', 'all', 'laneBarcode.html')
                    copyfile(demulti_stat_src, os.path.join(mfs_dest, 'laneBarcode.html'))
                    # Copy RunInfo.xml
                    run_info_xml_src = os.path.join(run.run_dir, 'RunInfo.xml')
                    if os.path.isfile(run_info_xml_src):
                        copyfile(run_info_xml_src, os.path.join(mfs_dest, 'RunInfo.xml'))
                    # Copy RunParameters.xml
                    run_parameters_xml_src = os.path.join(run.run_dir, 'RunParameters.xml')
                    if os.path.isfile(run_info_xml_src):
                        copyfile(run_parameters_xml_src, os.path.join(mfs_dest, 'RunParameters.xml'))
                    # Copy InterOp
                    interop_src = os.path.join(run.run_dir, 'InterOp')
                    if os.path.exists(interop_src):
                        copytree(interop_src, os.path.join(mfs_dest, 'InterOp'), dirs_exist_ok=True)
                except:
                    logger.warn('Could not copy demultiplex stats, InterOp metadata or XML files for run {}'.format(run.id))

            # Transfer to analysis server if flag is True
            if run.transfer_to_analysis_server:
                mail_recipients = CONFIG.get('mail', {}).get('recipients')
                logger.info('Transferring run {} to {} into {}'
                            .format(run.id,
                                    run.CONFIG['analysis_server']['host'],
                                    run.CONFIG['analysis_server']['sync']['data_archive']))
                run.transfer_run(transfer_file, mail_recipients)

            # Archive the run if indicated in the config file
            if 'storage' in CONFIG: #TODO: make sure archiving to PDC is not ongoing
                run.archive_run(CONFIG['storage']['archive_dirs'][run.sequencer_type])

    if run:
        # Determine the run type
        runObj = get_runObj(run, software)
        if not runObj:
            raise RuntimeError("Unrecognized instrument type or incorrect run folder {}".format(run))
        else:
            _process(runObj)
    else:
        data_dirs = CONFIG.get('analysis').get('data_dirs')
        for data_dir in data_dirs:
            # Run folder looks like DATE_*_*_*, the last section is the FC name.
            runs = glob.glob(os.path.join(data_dir, '[1-9]*_*_*_*'))
            for _run in runs:
                runObj = get_runObj(_run, software)
                if not runObj:
                    logger.warning('Unrecognized instrument type or incorrect run folder {}'.format(run))
                else:
                    try:
                        _process(runObj)
                    except:
                        # This function might throw and exception,
                        # it is better to continue processing other runs
                        logger.warning('There was an error processing the run {}'.format(run))
                        pass
