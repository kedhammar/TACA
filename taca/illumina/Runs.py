import os
import re
import csv
import logging
import subprocess
import shutil
import requests
from datetime import datetime

from taca.utils import misc
from taca.utils.misc import send_mail

from flowcell_parser.classes import RunParser

logger = logging.getLogger(__name__)

class Run(object):
    """ Defines an Illumina run
    """
    
    def __init__(self, run_dir, configuration):
        if not os.path.exists(run_dir):
            raise RuntimeError('Could not locate run directory {}'.format(run_dir))
        
        if 'analysis_server' not in configuration or \
            'bcl2fastq' not in configuration or \
            'samplesheets_dir' not in configuration:
            raise RuntimeError("configuration missing required entries "
                               "(analysis_server, bcl2fastq, samplesheets_dir)")
        
        if not os.path.exists(os.path.join(run_dir, 'runParameters.xml')) \
        and os.path.exists(os.path.join(run_dir, 'RunParameters.xml')):
            # In NextSeq runParameters is named RunParameters
            logger.warning("Renaming RunParameters.xml to runParameters.xml")
            os.rename(os.path.join(run_dir, 'RunParameters.xml'), os.path.join(run_dir, 'runParameters.xml'))
        elif not os.path.exists(os.path.join(run_dir, 'runParameters.xml')):
            raise RuntimeError('Could not locate runParameters.xml in run directory {}'.format(run_dir))
              
        self.run_dir = os.path.abspath(run_dir)
        self.id = os.path.basename(os.path.normpath(run_dir))
        pattern = r'(\d{6})_([ST-]*\w+\d+)_\d+_([AB]?)([A-Z0-9\-]+)'
        m = re.match(pattern, self.id)
        self.date = m.group(1)
        self.instrument = m.group(2)
        self.position = m.group(3)
        self.flowcell_id = m.group(4)
        self.CONFIG = configuration
        self._set_demux_folder(configuration)
        self._set_run_parser_obj(configuration) # get parser object to update DB
        # This flag tells TACA to move demultiplexed files to the analysis server
        self.transfer_to_analysis_server = True
        # Probably worth to add the samplesheet name as a variable too
        
    def demultiplex_run(self):
        raise NotImplementedError("Please Implement this method")


    def check_run_status(self):
        """
        This function checks the status of a run while in progress.
        In the case of HiSeq check that all demux have been done and in that case perform aggregation
        """
        run_dir    =  self.run_dir
        dex_status =  self.get_run_status()
        #in this case I have already finished all demux jobs and I have aggregate all stasts unded Demultiplexing
        if  dex_status == 'COMPLETED':
            return None
        #otherwise check the status of running demux
        #collect all samplesheets generated before
        samplesheets =  glob.glob(os.path.join(run_dir, "*_[0-9].csv")) # a single digit... this hipotesis should hold for a while
        allDemuxDone = True
        for samplesheet in samplesheets:
            #fetch the id of this demux job
            demux_id = os.path.splitext(os.path.split(samplesheet)[1])[0].split("_")[1]
            #demux folder is
            demux_folder = os.path.join(run_dir, "Demultiplexing_{}".format(demux_id))
            #check if this job is done
            if os.path.exists(os.path.join(run_dir, demux_folder, 'Stats', 'DemultiplexingStats.xml')):
                allDemuxDone = allDemuxDone and True
                logger.info("Sub-Demultiplexing in {} completed.".format(demux_folder))
            else:
                allDemuxDone = allDemuxDone and False
                logger.info("Sub-Demultiplexing in {} not completed yet.".format(demux_folder))
        #in this case, I need to aggreate in the Demultiplexing folder all the results
        if allDemuxDone:
            self._aggregate_demux_results()
            #now I can initialise the RunParser
            self.runParserObj = RunParser(self.run_dir)
            #and now I can rename undetermined if needed
            lanes = misc.return_unique([lanes['Lane'] for lanes in  self.runParserObj.samplesheet.data])
            samples_per_lane =  self.get_samples_per_lane()
            for lane in lanes:
                if self.is_unpooled_lane(lane):
                    self._rename_undet(lane, samples_per_lane)


    def _set_run_type(self):
        raise NotImplementedError("Please Implement this method")

    def get_run_type(self):
        if self.run_type:
            return self.run_type
        else:
            raise RuntimeError("run_type not yet available!!")

    def _set_sequencer_type(self, configuration):
        raise NotImplementedError("Please Implement this method")

    def _get_sequencer_type(self):
        if self.sequencer_type:
            return self.sequencer_type
        else:
            raise RuntimeError("sequencer_type not yet available!!")
    
    def _set_run_parser_obj(self, configuration):
        self.runParserObj = RunParser(self.run_dir)
        if self.runParserObj.obj:
            self.runParserObj.obj['DemultiplexConfig'] = \
            {'Setup': {'Software': configuration.get('bcl2fastq',{})}}

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
        """
            Locate and parse the samplesheet for a run. The idea is that there is a folder in
            samplesheet_folders that contains a samplesheet named flowecell_id.csv.
        """
        current_year = '20' + self.id[0:2]
        samplesheets_dir = os.path.join(self.CONFIG['samplesheets_dir'],
                                                current_year)
        ssname = os.path.join(samplesheets_dir, '{}.csv'.format(self.flowcell_id))
        if os.path.exists(ssname):
            return ssname
        else:
            raise RuntimeError("not able to find samplesheet {}.csv in {}".format(self.flowcell_id, self.CONFIG['samplesheets_dir']))


    def _is_demultiplexing_done(self):
        return os.path.exists(os.path.join(self.run_dir,
                                           self._get_demux_folder(), 
                                           'Stats',
                                           'DemultiplexingStats.xml'))

    def _is_demultiplexing_started(self):
        return os.path.exists(os.path.join(self.run_dir, self._get_demux_folder()))

    def _is_sequencing_done(self):
        return os.path.exists(os.path.join(self.run_dir, 'RTAComplete.txt'))

    def get_run_status(self):
        """ Return the status of the run, that is the trello card where it needs to be placed
        """
        demux_started = self._is_demultiplexing_started() # True if demux is ongoing
        demux_done = self._is_demultiplexing_done() # True if demux is done
        sequencing_done = self._is_sequencing_done() # True if sequencing is done
        if sequencing_done and demux_done:
            return 'COMPLETED' # run is done, transfer might be ongoing.
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
        # generate new ssparser (from the renamed smaplesheet)
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
            index2 = ""
            is_dual_index = False
            if data_entry.get('index'):
                index = data_entry['index']
                if index in "NoIndex": #special case for HiSeq when one sample is alone in a lane
                    index = ""
                is_dual_index = False # default for Xten
                if data_entry.get('index2'):
                    index2 = data_entry['index2']
                    is_dual_index = True
                #specific for HiSeq, will disapper once we will use bcl2fastq_2.17
                #index = data_entry['Index'].replace('-', '').replace('NoIndex', '')
            index_size  = len(index)
            index2_size = len(index2)
            # compute the basemask
            base_mask = self._compute_base_mask(runSetup, index_size, is_dual_index, index2_size)
            base_mask_string = "".join(base_mask)
            # prepare the dictionary
            if base_mask_string not in base_masks[lane]:
                # first time I find such base mask in this lane,
                base_masks[lane][base_mask_string] = {'base_mask':base_mask,
                                                      'data' : []}
            base_masks[lane][base_mask_string]['data'].append(data_entry)

        return base_masks

    def _compute_base_mask(self, runSetup, index_size, dual_index_sample, index2_size):
        """
            Assumptions:
                - if runSetup is of size 3, then single index run
                - if runSetup is of size 4, then dual index run
        """
        bm = []
        dual_index_run = False
        if len(runSetup) > 4:
            raise RuntimeError("when generating base_masks looks like there are"
                               " more than 4 reads in the RunSetup.xml")

        for read in runSetup:
            cycles = int(read['NumCycles'])
            if read['IsIndexedRead'] == 'N':
                bm.append('Y' + str(cycles))
            else:
                if index_size > cycles:
                    # the size of the index of the sample sheet is larger than the 
                    # one specified by RunInfo.xml, somethig must be wrong
                    raise RuntimeError("when generating base_masks found index in"  
                                       "samplesheet larger than the index specifed in RunInfo.xml")
                is_first_index_read = int(read['Number']) == 2
                # now prepare the base mask for this index read
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
                # when working on the second read index I need to know if the sample is dual index or not
                    if dual_index_sample:
                        i_remainder = cycles - index2_size
                        if i_remainder > 0:
                            bm.append('I' + str(index2_size) + 'N' + str(i_remainder))
                        else:
                            bm.append('I' + str(cycles))
                    else:
                    # if this sample is not dual index but the run is, 
                    # then I need to ignore the second index completely
                        bm.append('N' + str(cycles))
        return bm

    def transfer_run(self, t_file, analysis, mail_recipients=None):
        """ Transfer a run to the analysis server. Will add group R/W permissions to
            the run directory in the destination server so that the run can be processed
            by any user/account in that group (i.e a functional account...). 
            :param str t_file: File where to put the transfer information
            :param bool analysis: Trigger analysis on remote server
        """
        # TODO: check the run type and build the correct rsync command
        # The option -a implies -o and -g which is not the desired behaviour
        command_line = ['rsync', '-Lav', '--no-o', '--no-g']
        # Add R/W permissions to the group
        command_line.append('--chmod=g+rw')
        # This horrible thing here avoids data dup when we use multiple indexes in a lane/FC
        command_line.append("--exclude=Demultiplexing_*/*_*") 
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
            logger.error("Cannot create a file in {}. "
                         "Check the run name, and the permissions.".format(self.id))
            raise e
        started = ("Started transfer of run {} on {}".format(self.id, datetime.now()))
        logger.info(started)
        # In this particular case we want to capture the exception because we want
        # to delete the transfer file
        try:
           msge_text="I was about to transfer with this command \n{}".format(command_line)
           logger.info(msge_text)
           print msge_text
           # misc.call_external_command(command_line, with_log_files=True, 
           #                            prefix="", log_dir=self.run_dir)
        except subprocess.CalledProcessError as exception:
            os.remove(os.path.join(self.run_dir, 'transferring'))
            #Send an email notifying that the transfer failed
            runname = self.id
            sbt = ("Rsync of run {} failed".format(runname))
            msg= """ Rsync of data for run {run} has failed!
                Raised the following exception:     {e}
            """.format(run=runname, e=exception)
            if mail_recipients:
                send_mail(sbt, msg, mail_recipients)

            raise exception

        logger.info('Adding run {} to {}'.format(self.id, t_file))
        with open(t_file, 'a') as tranfer_file:
            tsv_writer = csv.writer(tranfer_file, delimiter='\t')
            tsv_writer.writerow([self.id, str(datetime.now())])
        os.remove(os.path.join(self.run_dir, 'transferring'))

        #Send an email notifying that the transfer was successful 
        runname = self.id
        sbt = ("Rsync of data for run {} to Irma has finished".format(runname))
        msg= """ Rsync of data for run {run} to Irma has finished!
                          
        The run is available at : https://genomics-status.scilifelab.se/flowcells/{run}
        """.format(run=runname)
        if mail_recipients:
            send_mail(sbt, msg, mail_recipients)


        if analysis:
            # This needs to pass the runtype (i.e., Xten or HiSeq) and start the correct pipeline
            self.trigger_analysis()
        
    def archive_run(self, destination):
        """ Move run to the archive folder
            :param str destination: the destination folder
        """
        if destination and os.path.isdir(destination):
            logger.info('archiving run {}'.format(self.id))
            shutil.move(self.run_dir, os.path.join(destination, self.id))
        else:
            logger.warning("Cannot move run to archive, destination does not exist")

    def trigger_analysis(self):
        """ Trigger the analysis of the flowcell in the analysis sever.
            :param str run_id: run/flowcell id
        """
        if not self.CONFIG.get('analysis_server', {}):
            logger.warn(("No configuration found for remote analysis server. "
                     "Not triggering analysis of {}"
                     .format(os.path.basename(self.id))))
        else:
            url = ("http://{host}:{port}/flowcell_analysis/{dir}"
                   .format(host=self.CONFIG['analysis_server']['host'],
                           port=self.CONFIG['analysis_server']['port'],
                           dir=os.path.basename(self.id)))
            params = {'path': self.CONFIG['analysis_server']['sync']['data_archive']}
            try:
                r = requests.get(url, params=params)
                if r.status_code != requests.status_codes.codes.OK:
                    logger.warn(("Something went wrong when triggering the "
                                 "analysis of {}. Please check the logfile "
                                 "and make sure to start the analysis!"
                                .format(os.path.basename(self.run_id))))
                else:
                    logger.info('Analysis of flowcell {} triggered in {}'
                                .format(os.path.basename(self.run_id),
                                self.CONFIG['analysis'][self.run_type]['analysis_server']['host']))
                    a_file = os.path.join(self.CONFIG['analysis'][self.run_type]['status_dir'], 'analysis.tsv')
                    with open(a_file, 'a') as analysis_file:
                        tsv_writer = csv.writer(analysis_file, delimiter='\t')
                        tsv_writer.writerow([os.path.basename(self.run_id), str(datetime.now())])
            except requests.exceptions.ConnectionError:
                logger.warn(("Something went wrong when triggering the analysis "
                            "of {}. Please check the logfile and make sure to "
                            "start the analysis!".format(os.path.basename(self.run_id))))

    def send_mail(self, msg, rcp):
        """ Sends mail about run completion
        """
        already_seen = False
        runname = self.id
        sj = "{}".format(runname)
        misc.send_mail(sj, msg, rcp)

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
                    # Rows have two columns: run and transfer date
                    if row[0] == os.path.basename(self.id):
                        return True
            if os.path.exists(os.path.join(self.run_dir, 'transferring')):
                return True
            return False
        except IOError:
            return False


    def is_unpooled_lane(self, lane):
        """
            :param lane: lane identifier
            :type lane: string
            :rtype: boolean
            :returns: True if the samplesheet has one entry for that lane, False otherwise
        """
        count = 0
        for l in self.runParserObj.samplesheet.data:
            if l['Lane'] == lane:
                count += 1
        return count == 1

    def is_unpooled_run(self):
        """
            :param ss: SampleSheet reader
            :type ss: flowcell_parser.XTenSampleSheet
            :rtype: boolean
            :returns: True if the samplesheet has one entry per lane, False otherwise
        """
        ar = []
        for l in self.runParserObj.samplesheet.data:
            ar.append(l['Lane'])
        return len(ar)==len(set(ar))

    def get_samples_per_lane(self):
        """
        :param ss: SampleSheet reader
        :type ss: flowcell_parser.XTenSampleSheet
        :rtype: dict
        :returns: dictionnary of lane:samplename
        """
        ss = self.runParserObj.samplesheet
        d={}
        for l in ss.data:
            s=l[ss.dfield_snm].replace("Sample_", "").replace("-", "_")
            d[l['Lane']]=l[ss.dfield_snm]

        return d



    def _rename_undet(self, lane, samples_per_lane):
        """Renames the Undetermined fastq file by prepending the sample name in front of it

        :param run: the path to the run folder
        :type run: str
        :param status: the demultiplexing status
        :type status: str
        :param samples_per_lane: lane:sample dict
        :type status: dict
        """
        run = self.run_dir
        dmux_folder = self.demux_dir
        for file in glob.glob(os.path.join(run, dmux_folder, "Undetermined*L0?{}*".format(lane))):
            old_name=os.path.basename(file)
            old_name_comps=old_name.split("_")
            old_name_comps[1]=old_name_comps[0]# replace S0 with Undetermined
            old_name_comps[0]=samples_per_lane[lane]#replace Undetermined with samplename
            for index, comp in enumerate(old_name_comps):
                if comp.startswith('L00'):
                    old_name_comps[index]=comp.replace('L00','L01')#adds a 1 as the second lane number in order to differentiate undetermined from normal in piper

            new_name="_".join(old_name_comps)
            logger.info("Renaming {} to {}".format(file, os.path.join(os.path.dirname(file), new_name)))
            os.rename(file, os.path.join(os.path.dirname(file), new_name))



