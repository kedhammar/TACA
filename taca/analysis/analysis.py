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
        if run.get_run_type() == 'NON-NGI-RUN':
            #For now MiSeq specific case. Process only NGI-run, skip all the others (PhD student runs)
            logger.warn('Run {} marked as NON-NGI-RUN, TACA will skip this and move the run to no-sync directory')
            run.archive_run(CONFIG['storage']['archive_dirs'][run.sequencer_type])
            return None
        
        if run.get_run_status() == 'SEQUENCING':
            # Check status files and say i.e Run in second read, maybe something
            # even more specific like cycle or something
            logger.info('Run {} is not finished yet'.format(run.id))
        elif  run.get_run_status() == 'TO_START':
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
        t_file = os.path.join(CONFIG['analysis']['status_dir'], 'transfer.tsv')
        if run.get_run_status() == 'COMPLETED' and not run.is_transferred(t_file):
            logger.info(("Preprocessing of run {} is finished, check if "
                             "run has been transferred and transfer it "
                             "otherwise".format(run.id)))
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
            _upload_to_statusdb(run)
            #if QC is ok tranfer the run in th appropriate server
            if run_QC_status or force_trasfer: #runs is not tranfer only it force_tranfer is False and QC_status false
                #tranfer the run, specify desitnation and if analysis needs to be started on the server
                logger.info("Run {} hasn't been transferred yet."
                                .format(run.id))
                logger.info('Transferring run {} to {} into {}'
                                    .format(run.id,
                                            run.CONFIG['analysis_server']['host'],
                                            run.CONFIG['analysis_server']['sync']['data_archive']))
                run.transfer_run(CONFIG['storage']['archive_dirs'][run.sequencer_type], t_file,  False) #do not trigger analysis
            else:
                logger.warn('Run {} failed qc, transferring will not take place'.format(run.id))
        elif run.is_transferred(t_file):
            logger.info('Run {} already transferred to analysis server, skipping it'.format(run.id))
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





