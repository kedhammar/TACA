import os
import re
import csv
import logging
import subprocess
import shutil
import glob
import json

from datetime import datetime

from taca.utils import misc
from taca.utils.misc import send_mail
from flowcell_parser.classes import RunParser, LaneBarcodeParser, SampleSheetParser

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
            logger.warning("Creating link from runParameters.xml to RunParameters.xml")
            os.symlink('RunParameters.xml', os.path.join(run_dir, 'runParameters.xml'))
        elif not os.path.exists(os.path.join(run_dir, 'runParameters.xml')):
            raise RuntimeError('Could not locate runParameters.xml in run directory {}'.format(run_dir))

        self.run_dir = os.path.abspath(run_dir)
        self.id = os.path.basename(os.path.normpath(run_dir))
        pattern = r'(\d{6,8})_([ST-]*\w+\d+)_\d+_([AB]?)([A-Z0-9\-]+)'
        m = re.match(pattern, self.id)
        self.date = m.group(1)
        self.instrument = m.group(2)
        self.position = m.group(3)
        self.flowcell_id = m.group(4)
        self.CONFIG = configuration
        self.demux_dir = "Demultiplexing"
        self.demux_summary = dict()
        self.runParserObj = RunParser(self.run_dir)
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
        # Check the status of running demux
        # Collect all samplesheets generated before
        samplesheets =  glob.glob(os.path.join(run_dir, "*_[0-9].csv")) # A single digit, this hypothesis should hold for a while
        all_demux_done = True
        for samplesheet in samplesheets:
            demux_id = os.path.splitext(os.path.split(samplesheet)[1])[0].split("_")[1]
            demux_folder = os.path.join(run_dir, "Demultiplexing_{}".format(demux_id))
            # Check if this job is done
            if os.path.exists(os.path.join(run_dir, demux_folder, 'Stats', 'DemultiplexingStats.xml')):
                all_demux_done = all_demux_done and True
                demux_log = os.path.join(run_dir, "demux_{}_bcl2fastq.err".format(demux_id))
                errors, warnings, error_and_warning_messages = self._check_demux_log(demux_id, demux_log)
                self.demux_summary[demux_id] = {'errors' : errors,
                                                'warnings' : warnings,
                                                'error_and_warning_messages' : error_and_warning_messages
                                               }
                if errors or warnings:
                    logger.info("Sub-Demultiplexing in {} completed with {} errors and {} warnings!".format(demux_folder, errors, warnings))
                else:
                    logger.info("Sub-Demultiplexing in {} completed without any error or warning.".format(demux_folder))
            else:
                all_demux_done = all_demux_done and False
                logger.info("Sub-Demultiplexing in {} not completed yet.".format(demux_folder))

        # All demux jobs finished and all stats aggregated under Demultiplexing
        # Aggreate all the results in the Demultiplexing folder
        if  all_demux_done and dex_status!='COMPLETED':
            dex_status = 'COMPLETED'
            self._aggregate_demux_results()
            self.runParserObj = RunParser(self.run_dir)
            # Rename undetermined if needed
            lanes = misc.return_unique([lanes['Lane'] for lanes in self.runParserObj.samplesheet.data])
            samples_per_lane = self.get_samples_per_lane()
            for lane in lanes:
                if self.is_unpooled_lane(lane):
                    self._rename_undet(lane, samples_per_lane)
            return None

    def _check_demux_log(self, demux_id, demux_log):
        """
        This function checks the log files of bcl2fastq
        Errors or warnings will be captured and email notifications will be sent
        """
        with open(demux_log, 'r') as demux_log_file:
            demux_log_content = demux_log_file.readlines()
            pattern = r'Processing completed with (\d+) errors and (\d+) warnings'
            match = re.search(pattern, demux_log_content[-1])
            if match:
                errors = int(match.group(1))
                warnings = int(match.group(2))
                error_and_warning_messages = []
                if errors or warnings:
                    for line in demux_log_content:
                        if 'ERROR' in line or 'WARN' in line:
                            error_and_warning_messages.append(line)
                return errors, warnings, error_and_warning_messages
            else:
                raise RuntimeError("Bad format with log file demux_{}_bcl2fastq.err".format(demux_id))

    def _set_run_type(self):
        raise NotImplementedError("Please Implement this method")

    def get_run_type(self):
        if self.run_type:
            return self.run_type
        else:
            raise RuntimeError("run_type not yet available!!")

    def _set_sequencer_type(self, configuration):
        raise NotImplementedError("Please Implement this method")

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
        try:
            # Only implemented for some, (e.g. NovaSeqXPlus)
            # Will raise AttributeError if not implemented.
            current_year = self._current_year()
        except AttributeError:
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
                                           'Stats.json'))

    def _is_demultiplexing_started(self):
        return os.path.exists(os.path.join(self.run_dir, self._get_demux_folder()))

    def _is_sequencing_done(self):
        return os.path.exists(os.path.join(self.run_dir, 'RTAComplete.txt')) and os.path.exists(os.path.join(self.run_dir, 'CopyComplete.txt'))

    def get_run_status(self):
        """ Return the current status of the run.
        """
        demux_started = self._is_demultiplexing_started()
        demux_done = self._is_demultiplexing_done()
        sequencing_done = self._is_sequencing_done()
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
                if index in ["NoIndex", "NOINDEX"]: #special case for HiSeq when one sample is alone in a lane
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
                            if index2_size == 0:
                                bm.append('N' + str(cycles)) #possible if same lane has single and dual index samples
                            else:
                                bm.append('I' + str(index2_size) + 'N' + str(i_remainder))
                        else:
                            bm.append('I' + str(cycles))
                    else:
                    # if this sample is not dual index but the run is,
                    # then I need to ignore the second index completely
                        bm.append('N' + str(cycles))
        return bm

    def transfer_run(self, t_file, mail_recipients=None):
        """ Transfer a run to the analysis server. Will add group R/W permissions to
            the run directory in the destination server so that the run can be processed
            by any user/account in that group (i.e a functional account...).
            :param str t_file: File where to put the transfer information
        """
        # The option -a implies -o and -g which is not the desired behaviour
        command_line = ['rsync', '-LtDrv']
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
           msge_text="I am about to transfer with this command \n{}".format(command_line)
           logger.info(msge_text)
           misc.call_external_command(command_line, with_log_files=True,
                                       prefix="", log_dir=self.run_dir)
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
        sbt = ("Rsync of data for run {} to the analysis cluster has finished".format(runname))
        msg= """ Rsync of data for run {run} to the analysis cluster has finished!

        The run is available at : https://genomics-status.scilifelab.se/flowcells/{run}
        """.format(run=runname)
        if mail_recipients:
            send_mail(sbt, msg, mail_recipients)

    def archive_run(self, destination):
        """ Move run to the archive folder
            :param str destination: the destination folder
        """
        if destination and os.path.isdir(destination):
            logger.info('archiving run {}'.format(self.id))
            shutil.move(self.run_dir, os.path.join(destination, self.id))
        else:
            logger.warning("Cannot move run to archive, destination does not exist")

    def send_mail(self, sbt, msg, rcp):
        """ Sends mail about run completion
        """
        already_seen = False
        runname = self.id
        if not sbt:
            sbt = "{}".format(runname)
        misc.send_mail(sbt, msg, rcp)

    def is_transferred(self, transfer_file):
        """ Checks wether a run has been transferred to the analysis server or not.
            Returns true in the case in which the tranfer is finished or ongoing.
            :param str transfer_file: Path to file with information about transferred runs
        """
        try:
            with open(transfer_file, 'r') as file_handle:
                transfer_file_contents = csv.reader(file_handle, delimiter='\t')
                for row in transfer_file_contents:
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

    def get_samples_per_lane(self):
        """
        :param ss: SampleSheet reader
        :type ss: flowcell_parser.XTenSampleSheet
        :rtype: dict
        :returns: dictionnary of lane:samplename
        """
        ss = self.runParserObj.samplesheet
        d = {}
        for l in ss.data:
            d[l['Lane']] = l[ss.dfield_snm]
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

    def _aggregate_demux_results_simple_complex(self):
        run_dir = self.run_dir
        runSetup = self.runParserObj.runinfo.get_read_configuration()
        demux_folder = os.path.join(self.run_dir , self.demux_dir)
        samplesheets = glob.glob(os.path.join(run_dir, "*_[0-9].csv"))

        index_cycles = [0, 0]
        for read in runSetup:
            if read['IsIndexedRead'] == 'Y':
                if int(read['Number']) == 2:
                    index_cycles[0] = int(read['NumCycles'])
                else:
                    index_cycles[1] = int(read['NumCycles'])

        # Prepare a list for lanes with NoIndex samples
        noindex_lanes = []
        for entry in self.runParserObj.samplesheet.data:
            if entry['index'].upper() == 'NOINDEX' or (entry['index'] == '' and entry['index2'] == ''):
                noindex_lanes.append(entry['Lane'])

        # Prepare a dict with the lane, demux_id and index_length info based on the sub-samplesheets
        # This is for the purpose of deciding simple_lanes and complex_lanes, plus we should start with the Stats.json file from which demux_id for each lane
        lane_demuxid_indexlength = dict()
        for samplesheet in samplesheets:
            demux_id = os.path.splitext(os.path.split(samplesheet)[1])[0].split("_")[1]
            ssparser = SampleSheetParser(samplesheet)
            for row in ssparser.data:
                if row['Lane'] not in lane_demuxid_indexlength.keys():
                    lane_demuxid_indexlength[row['Lane']] = {demux_id: [len(row.get('index','')), len(row.get('index2',''))]}
                elif demux_id not in lane_demuxid_indexlength[row['Lane']].keys():
                    lane_demuxid_indexlength[row['Lane']][demux_id] = [len(row.get('index','')), len(row.get('index2',''))]
                else:
                    pass

        simple_lanes = dict()
        complex_lanes = dict()
        for key, value in lane_demuxid_indexlength.items():
            if len(value) == 1:
                simple_lanes[key] = list(value.keys())[0]
            else:
                for vk, vv in value.items():
                    if key not in complex_lanes.keys():
                        complex_lanes[key] = {vk: vv}
                    else:
                        # Dual and longer indexes have higher priority
                        if 0 in list(complex_lanes[key].values())[0] and 0 not in vv:
                            complex_lanes[key] = {vk: vv}
                        elif (0 in list(complex_lanes[key].values())[0] and 0 in vv) or (0 not in list(complex_lanes[key].values())[0] and 0 not in vv):
                            if sum(vv) > sum(list(complex_lanes[key].values())[0]):
                                complex_lanes[key] = {vk: vv}
                        else:
                            pass

        # Case with only one sub-demultiplexing
        if len(complex_lanes) == 0 and len(samplesheets) == 1:
            demux_id_folder_name = "Demultiplexing_0" # in this case this is the only demux dir
            demux_id_folder = os.path.join(run_dir, demux_id_folder_name)
            # Special case that when we assign fake indexes for NoIndex samples
            if noindex_lanes and index_cycles != [0, 0]:
                # We first softlink the FastQ files of undet as the FastQ files of samples
                sample_counter = 1
                for entry in sorted(self.runParserObj.samplesheet.data, key=lambda k: k['Lane']):
                    lane = entry['Lane']
                    project = entry['Sample_Project']
                    sample = entry['Sample_ID']
                    project_dest = os.path.join(demux_folder, project)
                    if not os.path.exists(project_dest):
                        os.makedirs(project_dest)
                    sample_dest = os.path.join(project_dest, sample)
                    if not os.path.exists(sample_dest):
                        os.makedirs(sample_dest)
                    for file in glob.glob(os.path.join(demux_id_folder, "Undetermined*L0?{}*".format(lane))):
                        old_name = os.path.basename(file)
                        old_name_comps = old_name.split("_")
                        new_name_comps = [sample.replace('Sample_',''), 'S{}'.format(str(sample_counter))] + old_name_comps[2:]
                        new_name = "_".join(new_name_comps)
                        os.symlink(file, os.path.join(sample_dest, new_name))
                        logger.info("For undet sample {}, renaming {} to {}".format(sample.replace('Sample_',''), old_name, new_name))
                    sample_counter += 1
                # Make a softlink of lane.html
                html_report_lane_source = os.path.join(run_dir, demux_id_folder_name, "Reports", "html", self.flowcell_id, "all", "all", "all", "lane.html")
                html_report_lane_dest = os.path.join(demux_folder, "Reports", "html", self.flowcell_id, "all", "all", "all", "lane.html")
                if not os.path.isdir(os.path.dirname(html_report_lane_dest)):
                    os.makedirs(os.path.dirname(html_report_lane_dest))
                os.symlink(html_report_lane_source, html_report_lane_dest)

                # Modify the laneBarcode.html file
                html_report_laneBarcode = os.path.join(run_dir,
                                                       demux_id_folder_name,
                                                       "Reports",
                                                       "html",
                                                       self.flowcell_id,
                                                       "all",
                                                       "all",
                                                       "all",
                                                       "laneBarcode.html"
                                                       )
                html_report_laneBarcode_parser = LaneBarcodeParser(html_report_laneBarcode)
                lane_project_sample = dict()
                for entry in html_report_laneBarcode_parser.sample_data:
                    if entry['Sample'] != 'Undetermined':
                        lane_project_sample[entry['Lane']] = {'Project': entry['Project'],
                                                              'Sample': entry['Sample']
                                                              }
                for entry in html_report_laneBarcode_parser.sample_data[:]:
                    if entry['Sample'] == 'Undetermined':
                        entry['Project'] = lane_project_sample[entry['Lane']]['Project']
                        entry['Sample'] = lane_project_sample[entry['Lane']]['Sample']
                    else:
                        html_report_laneBarcode_parser.sample_data.remove(entry)
                html_report_laneBarcode_parser.sample_data = sorted(html_report_laneBarcode_parser.sample_data,
                                                                    key=lambda k: (k['Lane'].lower(), k['Sample']))
                new_html_report_laneBarcode = os.path.join(demux_folder,
                                                           "Reports",
                                                           "html",
                                                           self.flowcell_id,
                                                           "all",
                                                           "all",
                                                           "all",
                                                           "laneBarcode.html"
                                                           )
                _generate_lane_html(new_html_report_laneBarcode, html_report_laneBarcode_parser)

                if not os.path.exists(os.path.join(demux_folder, "Stats")):
                    os.makedirs(os.path.join(demux_folder, "Stats"))
                # Modify the Stats.json file
                stat_json_source = os.path.join(run_dir, demux_id_folder_name, "Stats", "Stats.json")
                stat_json_new = os.path.join(demux_folder, "Stats", "Stats.json")
                with open(stat_json_source) as json_data:
                    data = json.load(json_data)
                # Fix the sample stats per lane
                for entry in data['ConversionResults'][:]:
                    del entry['DemuxResults'][0]['IndexMetrics']
                    entry['DemuxResults'][0].update(entry['Undetermined'])
                    del entry['Undetermined']
                # Reset unknown barcodes list
                for entry in data['UnknownBarcodes'][:]:
                    entry['Barcodes'] = {'unknown': 1}
                # Write to a new Stats.json file
                with open(stat_json_new, 'w') as stat_json_new_file:
                    json.dump(data, stat_json_new_file)
            # This is the simple case, Demultiplexing dir is simply a symlink to the only sub-demultiplexing dir
            else:
                elements = [element for element  in  os.listdir(demux_id_folder) ]
                for element in elements:
                    if "Stats" not in element: #skip this folder and treat it differently to take into account the NoIndex case
                        source = os.path.join(demux_id_folder, element)
                        dest = os.path.join(self.run_dir, self.demux_dir, element)
                        os.symlink(source, dest)
                os.makedirs(os.path.join(self.run_dir, "Demultiplexing", "Stats"))
                # Fetch the lanes that have NoIndex
                statsFiles = glob.glob(os.path.join(demux_id_folder, "Stats", "*" ))
                for source in statsFiles:
                    source_name = os.path.split(source)[1]
                    if source_name not in ["DemultiplexingStats.xml", "AdapterTrimming.txt", "ConversionStats.xml", "Stats.json"]:
                        lane = os.path.splitext(os.path.split(source)[1])[0][-1] #lane
                        if lane not in noindex_lanes:
                            dest = os.path.join(self.run_dir, self.demux_dir, "Stats", source_name)
                            os.symlink(source, dest)
                for file in ["DemultiplexingStats.xml", "AdapterTrimming.txt", "ConversionStats.xml", "Stats.json"]:
                    source = os.path.join(self.run_dir, demux_id_folder_name, "Stats", file)
                    dest = os.path.join(self.run_dir, "Demultiplexing", "Stats", file)
                    os.symlink(source, dest)
            return True

        # Case with multiple sub-demultiplexings
        html_reports_lane = []
        html_reports_laneBarcode = []
        stats_json = []
        for samplesheet in samplesheets:
            ssparser = SampleSheetParser(samplesheet)
            demux_id = os.path.splitext(os.path.split(samplesheet)[1])[0].split("_")[1]
            demux_id_folder = os.path.join(run_dir, "Demultiplexing_{}".format(demux_id))
            html_report_lane = os.path.join(run_dir,
                                            "Demultiplexing_{}".format(demux_id),
                                            "Reports",
                                            "html",
                                            self.flowcell_id,
                                            "all",
                                            "all",
                                            "all",
                                            "lane.html"
                                            )
            if os.path.exists(html_report_lane):
                html_reports_lane.append(html_report_lane)
            else:
                raise RuntimeError("Not able to find html report {}: possible cause is problem in demultiplexing".format(html_report_lane))

            html_report_laneBarcode = os.path.join(run_dir,
                                                   "Demultiplexing_{}".format(demux_id),
                                                   "Reports",
                                                   "html",
                                                   self.flowcell_id,
                                                   "all",
                                                   "all",
                                                   "all",
                                                   "laneBarcode.html"
                                                   )
            if os.path.exists(html_report_laneBarcode):
                html_reports_laneBarcode.append(html_report_laneBarcode)
            else:
                raise RuntimeError("Not able to find html report {}: possible cause is problem in demultiplexing".format(html_report_laneBarcode))

            stat_json = os.path.join(run_dir, "Demultiplexing_{}".format(demux_id), "Stats", "Stats.json")
            if os.path.exists(stat_json):
                stats_json.append(stat_json)
            else:
                raise RuntimeError("Not able to find Stats.json report {}: possible cause is problem in demultiplexing".format(stat_json))

            # Aggregate fastq
            lanes_samples = dict()
            for row in ssparser.data:
                if row['Lane'] not in lanes_samples.keys():
                    lanes_samples[row['Lane']] = [row['Sample_Name']]
                else:
                    lanes_samples[row['Lane']].append(row['Sample_Name'])

            # Special case that when we assign fake indexes for NoIndex samples
            if (set(list(lanes_samples.keys())) & set(noindex_lanes)) and index_cycles != [0, 0]:
                sample_counter = 1
                for entry in sorted(ssparser.data, key=lambda k: k['Lane']):
                    lane = entry['Lane']
                    project = entry['Sample_Project']
                    sample = entry['Sample_ID']
                    project_dest = os.path.join(demux_folder, project)
                    if not os.path.exists(project_dest):
                        os.makedirs(project_dest)
                    sample_dest = os.path.join(project_dest, sample)
                    if not os.path.exists(sample_dest):
                        os.makedirs(sample_dest)
                    for file in glob.glob(os.path.join(demux_id_folder, "Undetermined*L0?{}*".format(lane))):
                        old_name = os.path.basename(file)
                        old_name_comps = old_name.split("_")
                        new_name_comps = [sample.replace('Sample_', ''), 'S{}'.format(str(sample_counter))] + old_name_comps[2:]
                        new_name = "_".join(new_name_comps)
                        os.symlink(file, os.path.join(sample_dest, new_name))
                        logger.info("For undet sample {}, renaming {} to {}".format(sample.replace('Sample_', ''), old_name, new_name))
                    sample_counter += 1
            # Ordinary cases
            else:
                projects = [project for project in  os.listdir(demux_id_folder) if os.path.isdir(os.path.join(demux_id_folder, project))]
                for project in projects:
                    if project in "Reports" or project in "Stats":
                        continue
                    project_source = os.path.join(demux_id_folder, project)
                    project_dest = os.path.join(demux_folder, project)
                    if not os.path.exists(project_dest):
                        # There might be project seqeunced with multiple index lengths
                        os.makedirs(project_dest)
                    samples = [sample for sample in os.listdir(project_source) if os.path.isdir(os.path.join(project_source, sample))]
                    for sample in samples:
                        sample_source = os.path.join(project_source, sample)
                        sample_dest = os.path.join(project_dest, sample)
                        if not os.path.exists(sample_dest):
                            # There should never be the same sample sequenced with different index length,
                            # however a sample might be pooled in several lanes and therefore sequenced using different samplesheets
                            os.makedirs(sample_dest)
                        fastqfiles =  glob.glob(os.path.join(sample_source, "*.fastq*"))
                        for fastqfile in fastqfiles:
                            os.symlink(fastqfile, os.path.join(sample_dest, os.path.split(fastqfile)[1]))
                # Copy fastq files for undetermined and the undetermined stats for simple lanes only
                lanes_in_sub_samplesheet = []
                header = ['[Header]','[Data]','FCID','Lane', 'Sample_ID', 'Sample_Name', 'Sample_Ref', 'index', 'index2', 'Description', 'Control', 'Recipe', 'Operator', 'Sample_Project']
                with open(samplesheet, mode='r') as sub_samplesheet_file:
                    sub_samplesheet_reader = csv.reader(sub_samplesheet_file)
                    for row in sub_samplesheet_reader:
                        if row[0] not in header:
                            lanes_in_sub_samplesheet.append(row[1])
                lanes_in_sub_samplesheet = list(set(lanes_in_sub_samplesheet))
                for lane in lanes_in_sub_samplesheet:
                    if lane in simple_lanes.keys():
                        undetermined_fastq_files = glob.glob(os.path.join(run_dir,
                                                                          "Demultiplexing_{}".format(demux_id),
                                                                          "Undetermined_S0_L00{}*.fastq*".format(lane))) # Contains only simple lanes undetermined
                        for fastqfile in undetermined_fastq_files:
                            os.symlink(fastqfile, os.path.join(demux_folder, os.path.split(fastqfile)[1]))
                        DemuxSummaryFiles = glob.glob(os.path.join(run_dir,
                                                                   "Demultiplexing_{}".format(demux_id),
                                                                   "Stats",
                                                                   "*L{}*txt".format(lane)))
                        if not os.path.exists(os.path.join(demux_folder, "Stats")):
                            os.makedirs(os.path.join(demux_folder, "Stats"))
                        for DemuxSummaryFile in DemuxSummaryFiles:
                            os.symlink(DemuxSummaryFile, os.path.join(demux_folder, "Stats", os.path.split(DemuxSummaryFile)[1]))

        # Create the html reports
        # Start with the lane
        html_report_lane_parser = None
        for next_html_report_lane in html_reports_lane:
            if html_report_lane_parser is None:
                html_report_lane_parser = LaneBarcodeParser(next_html_report_lane)
            else:
                lanesInReport = [Lane['Lane'] for Lane in html_report_lane_parser.sample_data]
                next_html_report_lane_parser = LaneBarcodeParser(next_html_report_lane)
                for entry in next_html_report_lane_parser.sample_data:
                    if not entry['Lane'] in lanesInReport:
                        # If this is a new lane not included before
                        html_report_lane_parser.sample_data.append(entry)
        # Now all lanes have been inserted

        # NumberReads for total lane cluster/yields and total sample cluster/yields
        NumberReads_Summary = dict()
        # The numbers in Flowcell Summary also need to be aggregated if multiple demultiplexing is done
        Clusters_Raw = 0
        Clusters_PF = 0
        Yield_Mbases = 0
        for entry in html_report_lane_parser.sample_data:
            # Update NumberReads for total lane clusters
            NumberReads_Summary[entry['Lane']] = {'total_lane_cluster': int(entry['PF Clusters'].replace(',', '')),
                                                  'total_lane_yield': int(entry['Yield (Mbases)'].replace(',', ''))}
            Clusters_Raw += int(int(entry['PF Clusters'].replace(',', '')) / float(entry['% PFClusters']) * 100)
            Clusters_PF += int(entry['PF Clusters'].replace(',', ''))
            Yield_Mbases += int(entry['Yield (Mbases)'].replace(',', ''))
            if entry['Lane'] in complex_lanes.keys():
                entry['% Perfectbarcode'] = None
                entry['% One mismatchbarcode'] = None
        # Update the values in Flowcell Summary
        html_report_lane_parser.flowcell_data['Clusters (Raw)'] = '{:,}'.format(Clusters_Raw)
        html_report_lane_parser.flowcell_data['Clusters(PF)'] = '{:,}'.format(Clusters_PF)
        html_report_lane_parser.flowcell_data['Yield (MBases)'] = '{:,}'.format(Yield_Mbases)
        # Add lanes not present in this demux
        # Create the new lane.html
        new_html_report_lane_dir = _create_folder_structure(demux_folder, ['Reports', 'html', self.flowcell_id, 'all', 'all', 'all'])
        new_html_report_lane = os.path.join(new_html_report_lane_dir, 'lane.html')
        _generate_lane_html(new_html_report_lane, html_report_lane_parser)

        # Generate the laneBarcode
        html_report_laneBarcode_parser = None
        for next_html_report_laneBarcode in html_reports_laneBarcode:
            if html_report_laneBarcode_parser is None:
                html_report_laneBarcode_parser = LaneBarcodeParser(next_html_report_laneBarcode)
            else:
                # No need to check samples occuring in more than one file as it would be spotted while softlinking
                next_html_report_laneBarcode_parser = LaneBarcodeParser(next_html_report_laneBarcode)
                for entry in next_html_report_laneBarcode_parser.sample_data:
                    html_report_laneBarcode_parser.sample_data.append(entry)
        # For complex lanes, set all numbers of undetermined to 0. And only keep one such entry
        constant_keys = ['Lane', 'Barcode sequence', 'Project', 'Sample']
        modified_complex_lanes = []
        for entry in html_report_laneBarcode_parser.sample_data:
            if entry['Lane'] in list(complex_lanes.keys()) and entry['Project'] in 'default':
                if entry['Lane'] not in modified_complex_lanes:
                    for key in entry.keys():
                        if key not in constant_keys:
                            entry[key] = '0'
                    modified_complex_lanes.append(entry['Lane'])
                else:
                    html_report_laneBarcode_parser.sample_data.remove(entry)

        # Update NumberReads for total sample yields
        for entry in html_report_laneBarcode_parser.sample_data:
            if 'total_sample_cluster' not in NumberReads_Summary[entry['Lane']].keys():
                NumberReads_Summary[entry['Lane']]['total_sample_cluster'] = 0
                NumberReads_Summary[entry['Lane']]['total_sample_yield'] = 0
                if entry['Project'] != 'default':
                    NumberReads_Summary[entry['Lane']]['total_sample_cluster'] += int(entry['PF Clusters'].replace(',', ''))
                    NumberReads_Summary[entry['Lane']]['total_sample_yield'] += int(entry['Yield (Mbases)'].replace(',', ''))
            else:
                if entry['Project'] != 'default':
                    NumberReads_Summary[entry['Lane']]['total_sample_cluster'] += int(entry['PF Clusters'].replace(',', ''))
                    NumberReads_Summary[entry['Lane']]['total_sample_yield'] += int(entry['Yield (Mbases)'].replace(',', ''))

        # Calculate the numbers clusters/yields of undet reads
        for key, value in NumberReads_Summary.items():
            value['undet_cluster'] = value['total_lane_cluster'] - value['total_sample_cluster']
            value['undet_yield'] = value['total_lane_yield'] - value['total_sample_yield']

        # Update the cluster/yield info of undet for complex lanes
        for entry in html_report_laneBarcode_parser.sample_data:
            if entry['Project'] == 'default' and entry['Lane'] in complex_lanes.keys():
                entry['PF Clusters'] = '{:,}'.format(NumberReads_Summary[entry['Lane']]['undet_cluster'])
                entry['Yield (Mbases)'] = '{:,}'.format(NumberReads_Summary[entry['Lane']]['undet_yield'])

        # Fix special case that when we assign fake indexes for NoIndex samples
        if noindex_lanes and index_cycles != [0, 0]:
            lane_project_sample = dict()
            for entry in html_report_laneBarcode_parser.sample_data:
                if entry['Lane'] in noindex_lanes and entry['Sample'] != 'Undetermined':
                    lane_project_sample[entry['Lane']] = {'Project': entry['Project'],
                                                          'Sample': entry['Sample']}
            for entry in html_report_laneBarcode_parser.sample_data[:]:
                if entry['Lane'] in noindex_lanes and entry['Sample'] == 'Undetermined':
                    entry['Project'] = lane_project_sample[entry['Lane']]['Project']
                    entry['Sample'] = lane_project_sample[entry['Lane']]['Sample']
                elif entry['Lane'] in noindex_lanes and entry['Sample'] != 'Undetermined':
                    html_report_laneBarcode_parser.sample_data.remove(entry)

        # Sort sample_data: first by lane then by sample ID
        html_report_laneBarcode_parser.sample_data = sorted(html_report_laneBarcode_parser.sample_data,
                                                            key=lambda k: (k['Lane'].lower(), k['Sample']))

        # Update the values in Flowcell Summary
        html_report_laneBarcode_parser.flowcell_data['Clusters (Raw)'] = '{:,}'.format(Clusters_Raw)
        html_report_laneBarcode_parser.flowcell_data['Clusters(PF)'] = '{:,}'.format(Clusters_PF)
        html_report_laneBarcode_parser.flowcell_data['Yield (MBases)'] = '{:,}'.format(Yield_Mbases)
        # Generate the new report for laneBarcode.html
        new_html_report_laneBarcode = os.path.join(new_html_report_lane_dir, 'laneBarcode.html')
        _generate_lane_html(new_html_report_laneBarcode, html_report_laneBarcode_parser)

        # Create the DemultiplexingStats.xml (empty it is here only to say thay demux is done)
        DemultiplexingStats_xml_dir = _create_folder_structure(demux_folder, ['Stats'])
        # For creating DemuxSummary.txt files for complex lanes
        DemuxSummaryFiles_complex_lanes = dict()
        # Generate the Stats.json
        with open(os.path.join(DemultiplexingStats_xml_dir, 'Stats.json'), 'w') as json_data_cumulative:
            stats_list = {}
            for stat_json in stats_json:
                demux_id = re.findall('Demultiplexing_([0-9])', stat_json)[0]
                with open(stat_json) as json_data_partial:
                    data = json.load(json_data_partial)
                    if len(stats_list) == 0:
                        # First time I do this
                        stats_list['RunNumber'] = data['RunNumber']
                        stats_list['Flowcell'] = data['Flowcell']
                        stats_list['RunId'] = data['RunId']
                        stats_list['ConversionResults'] = data['ConversionResults']
                        stats_list['ReadInfosForLanes'] = data['ReadInfosForLanes']
                        stats_list['UnknownBarcodes'] = []
                    else:
                        # Update only the importat fields
                        lanes_present_in_stats_json = [entry['LaneNumber'] for entry in stats_list['ConversionResults']]
                        for ReadInfosForLanes_lane in data['ReadInfosForLanes']:
                            if ReadInfosForLanes_lane['LaneNumber'] not in lanes_present_in_stats_json:
                                stats_list['ReadInfosForLanes'].extend([ReadInfosForLanes_lane])
                        for ConversionResults_lane in data['ConversionResults']:
                            if ConversionResults_lane['LaneNumber'] in lanes_present_in_stats_json and str(ConversionResults_lane['LaneNumber']) in complex_lanes.keys():
                                # For complex lanes, we set all stats to 0, except for read number and yield which will use values from NumberReads_Summary
                                ConversionResults_lane['Undetermined']['NumberReads'] = NumberReads_Summary[str(ConversionResults_lane['LaneNumber'])]['undet_cluster']
                                ConversionResults_lane['Undetermined']['Yield'] = NumberReads_Summary[str(ConversionResults_lane['LaneNumber'])]['undet_yield']*1000000
                                ConversionResults_lane['Undetermined']['ReadMetrics'][0]['QualityScoreSum'] = 0
                                ConversionResults_lane['Undetermined']['ReadMetrics'][0]['TrimmedBases'] = 0
                                ConversionResults_lane['Undetermined']['ReadMetrics'][0]['Yield'] = 0
                                ConversionResults_lane['Undetermined']['ReadMetrics'][0]['YieldQ30'] = 0
                                if len([r for r in self.runParserObj.runinfo.data['Reads'] if r['IsIndexedRead'] == 'N']) == 2:
                                    ConversionResults_lane['Undetermined']['ReadMetrics'][1]['QualityScoreSum'] = 0
                                    ConversionResults_lane['Undetermined']['ReadMetrics'][1]['TrimmedBases'] = 0
                                    ConversionResults_lane['Undetermined']['ReadMetrics'][1]['Yield'] = 0
                                    ConversionResults_lane['Undetermined']['ReadMetrics'][1]['YieldQ30'] = 0
                                # Find the list containing info for this lane #TODO: can lane_to_update be removed?
                                lane_to_update = [entry for entry in stats_list['ConversionResults'] if entry['LaneNumber'] == ConversionResults_lane['LaneNumber']][0]
                                lane_to_update['DemuxResults'].extend(ConversionResults_lane['DemuxResults'])
                                lane_to_update['Undetermined'] = ConversionResults_lane['Undetermined']
                            else:
                                stats_list['ConversionResults'].extend([ConversionResults_lane])

                    for unknown_barcode_lane in data['UnknownBarcodes']:
                        if str(unknown_barcode_lane['Lane']) in simple_lanes.keys():
                            stats_list['UnknownBarcodes'].extend([unknown_barcode_lane])
                        elif str(unknown_barcode_lane['Lane']) in complex_lanes.keys():
                            if list(complex_lanes[str(unknown_barcode_lane['Lane'])].keys())[0] == demux_id:
                                # First have the list of unknown indexes from the top priority demux run
                                full_list_unknownbarcodes = unknown_barcode_lane
                                # Remove the samples involved in the other samplesheets
                                for samplesheet in samplesheets:
                                    demux_id_ss = os.path.splitext(os.path.split(samplesheet)[1])[0].split("_")[1]
                                    if demux_id_ss != demux_id:
                                        ssparser = SampleSheetParser(samplesheet)
                                        ssparser_data_lane = [row for row in ssparser.data if row['Lane'] == str(unknown_barcode_lane['Lane'])]
                                        for row in ssparser_data_lane:
                                            sample_idx1 = row.get('index','')
                                            sample_idx2 = row.get('index2','')
                                            idx_copy = tuple(full_list_unknownbarcodes['Barcodes'].keys())
                                            for idx in idx_copy:
                                                unknownbarcode_idx1 = idx.split('+')[0] if '+' in idx else idx
                                                unknownbarcode_idx2 = idx.split('+')[1] if '+' in idx else ''
                                                if sample_idx1 and sample_idx2:
                                                    comparepart_idx1 = sample_idx1 if len(sample_idx1) <= len(unknownbarcode_idx1) else sample_idx1[:len(unknownbarcode_idx1)]
                                                    comparepart_idx2 = sample_idx2 if len(sample_idx2) <= len(unknownbarcode_idx2) else sample_idx2[:len(unknownbarcode_idx2)]
                                                    if comparepart_idx1 == unknownbarcode_idx1[:len(comparepart_idx1)] and comparepart_idx2 == unknownbarcode_idx2[:len(comparepart_idx2)]:
                                                        del full_list_unknownbarcodes['Barcodes'][idx]
                                                elif sample_idx1 and not sample_idx2:
                                                    comparepart_idx1 = sample_idx1 if len(sample_idx1) <= len(unknownbarcode_idx1) else sample_idx1[:len(unknownbarcode_idx1)]
                                                    if comparepart_idx1 == unknownbarcode_idx1[:len(comparepart_idx1)]:
                                                        del full_list_unknownbarcodes['Barcodes'][idx]
                                                elif not sample_idx1 and sample_idx2:
                                                    comparepart_idx2 = sample_idx2 if len(sample_idx2) <= len(unknownbarcode_idx1) else sample_idx2[:len(unknownbarcode_idx1)]
                                                    if comparepart_idx1 == unknownbarcode_idx1[:len(comparepart_idx2)]:
                                                        del full_list_unknownbarcodes['Barcodes'][idx]
                                stats_list['UnknownBarcodes'].extend([full_list_unknownbarcodes])
                                DemuxSummaryFiles_complex_lanes[str(unknown_barcode_lane['Lane'])] = full_list_unknownbarcodes
                            else:
                                pass

            # Fix special case that when we assign fake indexes for NoIndex samples
            if noindex_lanes and index_cycles != [0, 0]:
                for entry in stats_list['ConversionResults'][:]:
                    if str(entry['LaneNumber']) in noindex_lanes:
                        del entry['DemuxResults'][0]['IndexMetrics']
                        entry['DemuxResults'][0].update(entry['Undetermined'])
                        del entry['Undetermined']
                # Reset unknown barcodes list
                for entry in stats_list['UnknownBarcodes'][:]:
                    if str(entry['Lane']) in noindex_lanes:
                        entry['Barcodes'] = {'unknown': 1}

            # Write the final version of Stats.json file
            json.dump(stats_list, json_data_cumulative)

        # Create DemuxSummary.txt files for complex lanes
        if len(DemuxSummaryFiles_complex_lanes) > 0:
            for key, value in DemuxSummaryFiles_complex_lanes.items():
                with open(os.path.join(DemultiplexingStats_xml_dir, 'DemuxSummaryF1L{}.txt'.format(key)), 'w') as DemuxSummaryFile:
                    DemuxSummaryFile.write('### Most Popular Unknown Index Sequences\n')
                    DemuxSummaryFile.write('### Columns: Index_Sequence Hit_Count\n')
                    for idx, count in value['Barcodes'].items():
                        DemuxSummaryFile.write('{}\t{}\n'.format(idx, count))

        # Now the run is formally COMPLETED
        open(os.path.join(DemultiplexingStats_xml_dir, 'DemultiplexingStats.xml'), 'a').close()
        return True


def _create_folder_structure(root, dirs):
    """Creates a fodler stucture rooted in root usinf all dirs listed in dirs (a list)
    returns the path to the deepest directory
    """
    path = root
    for dir in dirs:
        path = os.path.join(path, dir)
        if not os.path.exists(path):
            os.makedirs(path)
    return path

def _generate_lane_html(html_file, html_report_lane_parser):
    with open(html_file, 'w') as html:
        # HEADER
        html.write('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">\n')
        html.write('<html xmlns:bcl2fastq>\n')
        html.write('<link rel="stylesheet" href="../../../../Report.css" type="text/css">\n')
        html.write('<body>\n')
        html.write('<table width="100%"><tr>\n')
        html.write('<td><p><p>C6L1WANXX /\n')
        html.write('        [all projects] /\n')
        html.write('        [all samples] /\n')
        html.write('        [all barcodes]</p></p></td>\n')
        html.write('<td><p align="right"><a href="../../../../FAKE/all/all/all/laneBarcode.html">show barcodes</a></p></td>\n')
        html.write('</tr></table>\n')
        # FLOWCELL SUMMARY TABLE
        html.write('<h2>Flowcell Summary</h2>\n')
        html.write('<table border="1" ID="ReportTable">\n')
        html.write('<tr>\n')
        fc_keys = sorted(list(html_report_lane_parser.flowcell_data.keys()))
        for key in fc_keys:
            html.write('<th>{}</th>\n'.format(key))
        html.write('</tr>\n')
        html.write('<tr>\n')
        for key in fc_keys:
            html.write('<td>{}</td>\n'.format(html_report_lane_parser.flowcell_data[key]))
        html.write('</tr>\n')
        html.write('</table>\n')
        # LANE SUMMARY TABLE
        html.write('<h2>Lane Summary</h2>\n')
        html.write('<table border="1" ID="ReportTable">\n')
        html.write('<tr>\n')
        lane_keys = sorted(list(html_report_lane_parser.sample_data[0].keys()))
        for key in lane_keys:
            html.write('<th>{}</th>\n'.format(key))
        html.write('</tr>\n')

        for sample in html_report_lane_parser.sample_data:
            html.write('<tr>\n')
            for key in lane_keys:
                html.write('<td>{}</td>\n'.format(sample[key]))
            html.write('</tr>\n')
        html.write('</table>\n')
        # FOOTER
        html.write('<p></p>\n')
        html.write('</body>\n')
        html.write('</html>\n')
