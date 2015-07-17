import os
import re
import csv
import glob
import datetime
import platform
import logging
import subprocess
import shutil

from datetime import datetime
from taca.utils import misc

from flowcell_parser.classes import RunParametersParser, SampleSheetParser, RunParser

logger = logging.getLogger(__name__)

class Run(object):
    """ Defines an Illumina run
    """

    def __init__(self, run_dir, configuration):
        if not os.path.exists(run_dir) or not \
                os.path.exists(os.path.join(run_dir, 'runParameters.xml')):
            raise RuntimeError('Could not locate run directory {}'.format(run_dir))
        if 'analysis_server' not in configuration or \
            'bcl2fastq' not in configuration or \
                'samplesheets_dir' not in configuration:
            raise RuntimeError('configuration missing required entries (analysis_server, bcl2fastq, samplesheets_dir')
        self.run_dir = os.path.abspath(run_dir)
        self.id = os.path.basename(os.path.normpath(run_dir))
        pattern = r'(\d{6})_([ST-]*\w+\d+)_\d+_([AB]?)([A-Z0-9\-]+)'
        m = re.match(pattern, self.id)
        self.date        = m.group(1)
        self.instrument  = m.group(2)
        self.position    = m.group(3)
        self.flowcell_id = m.group(4)
        self.CONFIG      = configuration
        self._set_demux_folder(configuration)
        self.runParserObj = RunParser(self.run_dir)


    def demultiplex_run():
        raise NotImplementedError("Please Implement this method")
    
    def check_run_status():
        raise NotImplementedError("Please Implement this method")
    
    def demux_done(self):
        raise NotImplementedError("Please Implement this method")
    
    def post_demux(self):
        raise NotImplementedError("Please Implement this method")

    
    def check_QC(self):
        raise NotImplementedError("Please Implement this method")

    

    def _set_sequencer_type(self, configuration):
        raise NotImplementedError("Please Implement this method")
    
    
    def _get_sequencer_type(self):
        if self.sequencer_type:
            return self.sequencer_type
        else:
            raise RuntimeError("sequencer_type not yet available!!")


    def _set_demux_folder(self, configuration):
        self.demux_dir = "Demultiplexing"
        for option in self.CONFIG['bcl2fastq']['options']:
            if isinstance(option, dict) and option.get('output-dir'):
                _demux_dir = option.get('output-dir')


    def _get_demux_folder(self):
        if self.demux_dir:
            return self.demux_dir
        else:
            raise RuntimeError("demux_folder not yet available!!")



    def _get_samplesheet(self):
        raise NotImplementedError("Please Implement this method")



    def _is_demultiplexing_done(self):
        if os.path.exists(os.path.join(self.run_dir, self._get_demux_folder(), 'Stats', 'DemultiplexingStats.xml')):
            return True
        else:
            return False



    def _is_demultiplexing_started(self):
        if os.path.exists(os.path.join(self.run_dir, self._get_demux_folder())):
            return True
        else:
            return False


    def _is_sequencing_done(self):
        if os.path.exists(os.path.join(self.run_dir, 'RTAComplete.txt')):
            return True
        else:
            return False


    def get_run_status(self):
        """
        return the status of the run, that is the trello card where it needs to be placed
        """
        demux_started   = self._is_demultiplexing_started() # True if demux is ongoing
        demux_done      = self._is_demultiplexing_done()    # True if demux is done
        sequencing_done = self._is_sequencing_done()        # True if sequencing is done
        if sequencing_done and demux_done:
            return 'COMPLETED' # run is done, tranfer might beongoing.
        elif sequencing_done and demux_started and not demux_done:
            return 'IN_PROGRESS'
        elif sequencing_done and not demux_started:
            return 'TO_START'
        elif not sequencing_done:
            return 'SEQUENCING'
        else:
            raise RuntimeError('Unexpected status in get_run_status')




    def _generate_per_lane_base_mask(self):
        """
        This functions generate the base masks for each lane. This is required as often we
        run FC with lanes that have different index lengths. A more extreme case are lanes where
        multiple index lengths are present in the same lane.
        Hypotesis:
            - RunInfo.xml contains the configuration
            - this object contains a properly parsed samplesheet
            
        It returns an dict with a key for each lane:
        {lane1:
            {base_mask_string (e.g., Y150I6N2N8Y150):
                [ base_mask , [SampleSheetEntries]]
            }
            {base_mask_string (e.g., Y150I6N2I6N2Y150):
                [ base_mask , [SampleSheetEntries]]
            }
         lane2:
        }
            
        """
        #generate new ssparser (from the renamed smaplesheet)
        runSetup = self.runParserObj.runinfo.get_read_configuration()
        base_masks = {}
        if not self.runParserObj.samplesheet:
            raise RuntimeError("samplesheet not yet initialised")
        
        for data_entry in self.runParserObj.samplesheet.data:
            ## for each data_entry in my samplesheet (i.e., for each sample)
            lane  = data_entry['Lane']
            if lane not in base_masks:
                base_masks[lane] = {}
            index = ""
            is_dual_index = False
            if data_entry.get('index'):
                index = data_entry['index']
                is_dual_index = "-" in  data_entry['index']
            elif data_entry.get('Index'):
                #specific for HiSeq, will disapper once we will use bcl2fastq_2.17
                index = data_entry['Index'].replace('-', '').replace('NoIndex', '')
                is_dual_index = "-" in  data_entry['Index']
            index_size = len(index)
            ##I have computed the index size, and I know if I have a single or dual index
            if is_dual_index:
                index_size = index_size/2
            #compute the basemask
            base_mask = self._compute_base_mask(runSetup, index_size, is_dual_index)
            base_mask_string = "".join(base_mask)
            #prepare the dictionary
            if base_mask_string not in base_masks[lane]:
                #first time I find such base mask in this lane,
                base_masks[lane][base_mask_string] = {'base_mask':base_mask,
                                                      'data' : []}
            base_masks[lane][base_mask_string]['data'].append(data_entry)

        return base_masks




    def _compute_base_mask(self, runSetup, index_size, dual_index_sample):
        """
            Assumptions:
                - if runSetup is of size 3, then single index run
                - if runSetup is of size 4, then dual index run
        """
        bm = []
        dual_index_run = False
        if len(runSetup) != 3 and len(runSetup) != 4:
            raise RuntimeError("when generating base_masks flooks like there are more than 3 or 4 reads in the RunSetup.xml")

        for read in runSetup:
            cycles = int(read['NumCycles'])
            if read['IsIndexedRead'] == 'N':
                bm.append('Y' + str(cycles))
            else:
                if index_size > cycles:
                    #the size of the index of the sample sheet is larger than the one specified by RunInfo.xml, somethig must be wrong
                    raise RuntimeError("when generating base_masks found index in samplesheet larger than the index specifed in RunInfo.xml")
                is_first_index_read = int(read['Number']) == 2
                #now prepare the base mask for this index read
                if is_first_index_read:
                    i_remainder = cycles - index_size
                    if i_remainder > 0:
                        if index_size == 0:
                            bm.append('N' + str(cycles)) #special case (NoIndex)
                        else:
                            bm.append('I' + str(index_size) + 'N' + str(i_remainder))
                    else:
                        bm.append('I' + str(cycles))
                else:
                #when working on the second read index I need to know if the sample is dual index or not
                    if dual_index_sample:
                        i_remainder = cycles - index_size
                        if i_remainder > 0:
                            bm.append('I' + str(index_size) + 'N' + str(i_remainder))
                        else:
                            bm.append('I' + str(cycles))
                    else:
                    #if this sample is not dual index but the run is, then I need to ignore the second index completely
                            bm.append('N' + str(cycles))
        return bm



    def tranfer_run(self, destination, t_file, analysis):
        """ Transfer a run to the analysis server. Will add group R/W permissions to
            the run directory in the destination server so that the run can be processed
            by any user/account in that group (i.e a functional account...). Run will be
            moved to data_dir/nosync after transferred.
        
            :param str run: Run directory
            :param bool analysis: Trigger analysis on remote server
        """
        #TODO: chekc the run type and build the correct rsync command
        command_line = ['rsync', '-av']
        # Add R/W permissions to the group
        command_line.append('--chmod=g+rw')
        # rsync works in a really funny way, if you don't understand this, refer to
        # this note: http://silentorbit.com/notes/2013/08/rsync-by-extension/
        command_line.append("--include=*/")
        for to_include in self.CONFIG['analysis_server']['sync']['include']:
            command_line.append("--include={}".format(to_include))
        command_line.extend(["--exclude=*", "--prune-empty-dirs"])
        r_user = self.CONFIG['analysis_server']['user']
        r_host = self.CONFIG['analysis_server']['host']
        r_dir = self.CONFIG['analysis_server']['sync']['data_archive']
        remote = "{}@{}:{}".format(r_user, r_host, r_dir)
        command_line.extend([self.run_dir, remote])

        # Create temp file indicating that the run is being transferred
        try:
            open(os.path.join(self.run_dir, 'transferring'), 'w').close()
        except IOError as e:
            logger.error("Cannot create a file in {}. Check the run name, and the permissions.".format(self.id))
            raise e
        started = ("Started transfer of run {} on {}".format(self.id, datetime.now()))
        logger.info(started)
        # In this particular case we want to capture the exception because we want
        # to delete the transfer file
        try:
            misc.call_external_command(command_line, with_log_files=True, prefix="", log_dir=self.run_dir)
        except subprocess.CalledProcessError as exception:
            os.remove(os.path.join(self.run_dir, 'transferring'))
            raise exception

        logger.info('Adding run {} to {}'.format(self.id, t_file))
        with open(t_file, 'a') as tranfer_file:
            tsv_writer = csv.writer(tranfer_file, delimiter='\t')
            tsv_writer.writerow([run.id, str(datetime.now())])
        os.remove(os.path.join(self.run_dir, 'transferring'))

        #Now, let's move the run to nosync
        self.archive_run(destination)
    
        if analysis:
            #This needs to pass the runtype (i.e., Xten or HiSeq) and start the correct pipeline
            self.trigger_analysis()


    def archive_run(self, destination):
        if destination:
            logger.info('archiving run {}'.format(run))
            shutil.move(os.path.abspath(run), os.path.join(destination, os.path.basename(os.path.abspath(run))))


    def trigger_analysis(self):
        """ Trigger the analysis of the flowcell in the analysis sever.
        
            :param str run_id: run/flowcell id
        """
        if not CONFIG.get('analysis', {}).get(run_type,{}).get('analysis_server', {}):
            logger.warn(("No configuration found for remote analysis server. "
                     "Not triggering analysis of {}"
                     .format(os.path.basename(run_id))))
        else:
            url = ("http://{host}:{port}/flowcell_analysis/{dir}"
                        .format(host=CONFIG['analysis'][run_type]['analysis_server']['host'],
                       port=CONFIG['analysis'][run_type]['analysis_server']['port'],
                       dir=os.path.basename(run_id)))
            params = {'path': CONFIG['analysis'][run_type]['analysis_server']['sync']['data_archive']}
            try:
                r = requests.get(url, params=params)
                if r.status_code != requests.status_codes.codes.OK:
                    logger.warn(("Something went wrong when triggering the "
                                    "analysis of {}. Please check the logfile "
                                    "and make sure to start the analysis!"
                                    .format(os.path.basename(run_id))))
                else:
                    logger.info('Analysis of flowcell {} triggered in {}'
                                    .format(os.path.basename(run_id),
                                    CONFIG['analysis'][run_type]['analysis_server']['host']))
                    a_file = os.path.join(CONFIG['analysis'][run_type]['status_dir'], 'analysis.tsv')
                    with open(a_file, 'a') as analysis_file:
                        tsv_writer = csv.writer(analysis_file, delimiter='\t')
                        tsv_writer.writerow([os.path.basename(run_id), str(datetime.now())])
            except requests.exceptions.ConnectionError:
                logger.warn(("Something went wrong when triggering the analysis "
                            "of {}. Please check the logfile and make sure to "
                            "start the analysis!".format(os.path.basename(run_id))))





    def post_qc(self, qc_file, status):
        """ Checks wether a run has passed the final qc.
            :param str run: Run directory
            :param str qc_file: Path to file with information about transferred runs
        """
        already_seen=False
        runname=os.path.basename(os.path.abspath(self.run_dir))
        shortrun=runname.split('_')[0] + '_' +runname.split('_')[-1]
        with open(qc_file, 'ab+') as f:
            f.seek(0)
            for row in f:
                #Rows have two columns: run and transfer date
                if row.split('\t')[0] == runname:
                    already_seen=True
    
            if not already_seen:
                if status:
                    f.write("{}\tPASSED\n".format(runname))
                else:
                    sj="{} failed QC".format(runname)
                    cnt="""The run {run} has failed qc and will NOT be transfered to Nestor.
                    
                        The run might be available at : https://genomics-status.scilifelab.se/flowcells/{shortfc}
                    
                        To read the logs, run the following command on {server}
                        grep -A30 "Checking run {run}" {log}
                    
                        To force the transfer :
                        taca analysis transfer {rundir} """.format(run=runname, shortfc=shortrun, log=CONFIG['log']['file'], server=os.uname()[1], rundir=run)
                    rcp=CONFIG['mail']['recipients']
                    misc.send_mail(sj, cnt, rcp)
                    f.write("{}\tFAILED\n".format(os.path.basename(run)))


    def is_transferred(self, transfer_file):
        """ Checks wether a run has been transferred to the analysis server or not.
        Returns true in the case in which the tranfer is ongoing.
        
        :param str run: Run directory
        :param str transfer_file: Path to file with information about transferred runs
        """
        try:
            with open(transfer_file, 'r') as file_handle:
                t_f = csv.reader(file_handle, delimiter='\t')
                for row in t_f:
                    #Rows have two columns: run and transfer date
                    if row[0] == os.path.basename(self.run_id):
                        return True
            if os.path.exists(os.path.join(self.run_id, 'transferring')):
                return True
            return False
        except IOError:
            return False












