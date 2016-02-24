""" Analysis methods for TACA """
import csv
import glob
import logging
import os
import re


from taca.illumina.Runs import Run
from taca.illumina.HiSeqX_Runs import HiSeqX_Run
from taca.illumina.HiSeq_Runs import HiSeq_Run
from taca.illumina.MiSeq_Runs import MiSeq_Run
from taca.utils.config import CONFIG
from shutil import copyfile

import flowcell_parser.db as fcpdb
from   flowcell_parser.classes import RunParametersParser

logger = logging.getLogger(__name__)







def _run_type(run):
    """Tries to read runParameters.xml and returns the run type.
        
        :param run: run name identifier
        :type run: string
        :rtype: String
        :returns: returns a string with the sequencer type name, None if the sequencer type is unknown
    """
    
    rppath=os.path.join(run, 'runParameters.xml')
    try:
        rp=RunParametersParser(os.path.join(run, 'runParameters.xml'))
    except OSError:
        logger.warn("Cannot find the runParameters.xml file at {}. This is quite unexpected. please archive the run {} manually".format(rppath, run))
    else:
        #this information about the run type (with HiSeq2.5 applicationaName does not work anymore, but as for a long time we will have instruments not updated I need to find out somehting that works
        try:
            #Works for recent control software
            runtype=rp.data['RunParameters']["Setup"]["Flowcell"]
        except KeyError:
            #use this as second resource but print a warning in the logs
            logger.warn("Parsing runParameters to fecth instrument type, not found Flowcell information in it. Using ApplicaiotnName")
            runtype=rp.data['RunParameters']["Setup"].get("ApplicationName", "") # here makes sense to use get with default value "" -> so that it doesnt raise an exception in the next lines (in case ApplicationName is not found, get returns None)

        if "HiSeq X" in runtype:
            return 'HiSeqX'
        elif "MiSeq" in runtype:
            return 'MiSeq'
        elif "HiSeq" in runtype or "TruSeq" in runtype:
            return 'HiSeq'
        else:
            logger.warn("unrecognized runtype {}, cannot archive the run {}. Someone as likely bought a new sequencer without telling it to the bioinfo team".format(runtype, run))
            return None
    return None




def upload_to_statusdb(run_dir):
    """Function to upload run_dir informations to statusDB directly from click interface
        
        :param run_dir: run name identifier
        :type run: string
        :rtype: None
    """
    sequencer_type = _run_type(run_dir)
    if sequencer_type is 'HiSeqX':
        runObj = HiSeqX_Run(run_dir, CONFIG["analysis"]["HiSeqX"])
    elif sequencer_type is 'HiSeq':
        runObj = HiSeq_Run(run_dir, CONFIG["analysis"]["HiSeq"])
    elif sequencer_type is 'MiSeq':
        runObj = MiSeq_Run(run_dir, CONFIG["analysis"]["MiSeq"])
    _upload_to_statusdb(runObj)
    return None

def _upload_to_statusdb(run):
    """Triggers the upload to statusdb using the dependency flowcell_parser
    
    :param Run run: the object run
    """
    couch = fcpdb.setupServer(CONFIG)
    db    = couch[CONFIG['statusdb']['xten_db']]
    parser = run.runParserObj
    for element in parser.obj['samplesheet_csv']:
        if 'NoIndex' in element['index'] or not element['index']: #NoIndex in the case of HiSeq, empty in the case of HiSeqX
            lane = element['Lane'] # this is a lane with NoIndex
            #in this case PF Cluster is the number of undetermined reads, as they are all undet
            try:
                PFclusters = parser.obj['Undetermined'][lane]['unknown']
            except KeyError:
                logger.error("While taking extra care of lane {} of NoIndex type I found out that not all values were available".format(lane))
                continue
            #in Lanes_stats fix the lane yield
            parser.obj['illumina']['Demultiplex_Stats']['Lanes_stats'][int(lane) - 1]['PF Clusters'] = str(PFclusters)
            #now fix Barcode lane stats
            updated = 0 #check that only one update is made
            for sample in parser.obj['illumina']['Demultiplex_Stats']['Barcode_lane_statistics']:
                if lane in sample['Lane']:
                    updated +=1
                    sample['PF Clusters'] = str(PFclusters)
            if updated != 1:
                logger.error("While taking extra care of lane {} of NoIndex type I updated more than once the barcode_lane. This is too much to continue so I will fail.".format(lane))
                sys.exit()
            #If I am here it means I changed the html representation to somthing else to accomodate the wired things we do
            #someone told me that in such cases it is better to put a place holder for this
            parser.obj['illumina']['Demultiplex_Stats']['NotOriginal'] = "True"
    fcpdb.update_doc( db , parser.obj, over_write_db_entry=True)
    return None

def transfer_run(run_dir, analysis):
    """interface for click to tranfer a run to uppmax
    
    :param: string run_dir: the run to tranfer
    :param bool analysis: if trigger or not the analysis
    """
    sequencer_type = _run_type(run_dir)
    if sequencer_type is 'HiSeqX':
        runObj = HiSeqX_Run(run_dir, CONFIG["analysis"]["HiSeqX"])
    elif sequencer_type is 'HiSeq':
        runObj = HiSeq_Run(run_dir, CONFIG["analysis"]["HiSeq"])
    elif sequencer_type is 'MiSeq':
        runObj = MiSeq_Run(run_dir, CONFIG["analysis"]["MiSeq"])
    else:
        logger.error("looks like we bough a new sequencer and no-body told me about it...")
        return None
    runObj.transfer_run("nosync", os.path.join(CONFIG['analysis']['status_dir'], 'transfer.tsv'),
                       analysis) #do not start analsysis automatically if I force the tranfer
    return None




def run_preprocessing(run, force_trasfer=True):
    """Run demultiplexing in all data directories

    :param str run: Process a particular run instead of looking for runs
    :param bool force_tranfer: if set to True the FC is tranferred also if fails QC
    """
  

    def _process(run, force_trasfer):
        """Process a run/flowcell and transfer to analysis server

        :param taca.illumina.Run run: Run to be processed and transferred
        """
        if not run:
            #this is in case the methods are not yet implemented
            return None
        logger.info('Checking run {}'.format(run.id))
        
        t_file = os.path.join(CONFIG['analysis']['status_dir'], 'transfer.tsv')
        if run.is_transferred(t_file):
            #in this case I am either processing a run that is in tranfer or that has been already tranferred. Do nothing.
            #time to time this situation is due to runs that are copied back from NAS after a reboot. This check avoid faiulures
            logger.info('Run {} already transferred to analysis server, skipping it'.format(run.id))
            return None
        
        if run.get_run_status() == 'SEQUENCING':
            # Check status files and say i.e Run in second read, maybe something
            # even more specific like cycle or something
            logger.info('Run {} is not finished yet'.format(run.id))
        elif  run.get_run_status() == 'TO_START':
            if run.get_run_type() == 'NON-NGI-RUN':
                #For now MiSeq specific case. Process only NGI-run, skip all the others (PhD student runs)
                logger.warn('Run {} marked as {}, TACA will skip this and move the run to no-sync directory'.format(run.id, run.get_run_type()))
                run.archive_run(CONFIG['storage']['archive_dirs'][run.sequencer_type])
                return None
            #otherwhise it is fine, process it
            logger.info(("Starting BCL to FASTQ conversion and "
                             "demultiplexing for run {}".format(run.id)))
            run.demultiplex_run()

        elif run.get_run_status() == 'IN_PROGRESS':
            logger.info(("BCL conversion and demultiplexing process in "
                             "progress for run {}, skipping it"
                             .format(run.id)))
            #in the case of Xten retruns, in future have a look to Cycles.txt
            #in the case of HiSeq check that partial demux are done and performs aggregation if this is the case
            run.check_run_status()

        # previous elif might change the status to COMPLETED (in HiSeq), therefore to avoid skip
        # a cicle take the last if out of the elif
        if run.get_run_status() == 'COMPLETED':
            logger.info(("Preprocessing of run {} is finished, "
                             "tranferring it".format(run.id)))
            #in the case of of HiSeq this function computes undetermined indexes for NoIndex lanes
            if not run.compute_undetermined():
                return None
            #otherwise I can procced to QC
            #check the run QC
            run_QC_status = run.check_QC()
            #store QC results in appropriate file and mail user if failed
            qc_file = os.path.join(CONFIG['analysis']['status_dir'], 'qc.tsv')
            #this method is implemented in Runs
            run.post_qc(qc_file, run_QC_status, log_file=CONFIG['log']['file'], rcp=CONFIG['mail']['recipients'])
            #upload to statusDB
            if 'statusdb' in CONFIG:
                _upload_to_statusdb(run)
            #copy demultiplex stat file to shard file system for LIMS purpose
            try:
                mfs_dest = os.path.join(CONFIG['mfs_path'],"{}_data".format(_run_type(run).lower()),run.id, 'laneBarcode.html')
                logger.info('Copying demultiplex stat for run {} to {}'.format(run.id, CONFIG[]))
                if not os.path.exists(mfs_dest):
                    os.mkdir(mfs_dest)
                demulti_stat_src = os.path.join(run.run_dir, run.demux_dir, 'Reports', 'html', run.flowcell_id, 'all', 'all', 'all', 'laneBarcode.html')
                copyfile(demulti_stat_src, mfs_dest)
            except:
                logger.warn('Could not copy demultiplex stat file for run {}'.format(run.id))
            #tranfer data to uppmax
            logger.info('Transferring run {} to {} into {}'
                                .format(run.id,
                                run.CONFIG['analysis_server']['host'],
                                run.CONFIG['analysis_server']['sync']['data_archive']))
            run.transfer_run(CONFIG['storage']['archive_dirs'][run.sequencer_type], t_file,  False) #do not trigger analysis
            
        return None

    if run:
        #needs to guess what run type I have (HiSeq, MiSeq, HiSeqX)
        sequencer_type = _run_type(run)
        if sequencer_type is 'HiSeqX':
            runObj = HiSeqX_Run(run, CONFIG["analysis"]["HiSeqX"])
        elif sequencer_type is 'HiSeq':
            runObj = HiSeq_Run(run, CONFIG["analysis"]["HiSeq"])
        elif sequencer_type is 'MiSeq':
            runObj = MiSeq_Run(run, CONFIG["analysis"]["MiSeq"])
        else:
            raise RuntimeError("New instrument type {}".format(sequencer_type))
        _process(runObj, force_trasfer)
    else:
        data_dirs = CONFIG.get('analysis').get('data_dirs')
        for data_dir in data_dirs:
            runs = glob.glob(os.path.join(data_dir, '1*XX'))
            # Try MiSeq runs as well
            if not runs:
                runs = glob.glob(os.path.join(data_dir, '1*000000000*'))
            for _run in runs:
                sequencer_type = _run_type(_run)
                if sequencer_type is 'HiSeqX':
                    runObj = HiSeqX_Run(_run, CONFIG["analysis"]["HiSeqX"])
                elif sequencer_type is 'HiSeq':
                    runObj = HiSeq_Run(_run, CONFIG["analysis"]["HiSeq"])
                elif sequencer_type is 'MiSeq':
                    runObj = MiSeq_Run(_run, CONFIG["analysis"]["MiSeq"])
                else:
                    raise RuntimeError("New instrument type {}".format(sequencer_type))
                _process(runObj, force_trasfer)





