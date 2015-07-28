import os
import re
import csv
import glob
import shutil
from datetime import datetime
from taca.utils.filesystem import chdir, control_fastq_filename
from taca.illumina.Runs import Run
from taca.utils import misc
from flowcell_parser.classes import RunParametersParser, SampleSheetParser, RunParser, LaneBarcodeParser


import logging

logger = logging.getLogger(__name__)

class HiSeq_Run(Run):

    def __init__(self,  path_to_run, samplesheet_folders):
        super(HiSeq_Run, self).__init__( path_to_run, samplesheet_folders)

    def _sequencer_type(self):
            return "HiSeq"


    def get_run_info(self):
        """
        Parse the RunInfo.xml file into a dict
        """
        f = os.path.join(self.path,'RunInfo.xml')
        if not os.path.exists(f):
            return {}
        with open(f) as fh:
            rip = RunInfoParser()
            runinfo = rip.parse(fh)
        return runinfo


    def get_projects(self):
        ssheet = self.samplesheet
        if ssheet is None:
            return None
        samplesheet_Obj = HiSeqSampleSheet(ssheet)
        return samplesheet_Obj.return_projects()



    def get_run_mode(self):
        if self.runParserObj:
            if self.runParserObj.runparameters.data.has_key('RunParameters') and \
               self.runParserObj.runparameters.data['RunParameters'].has_key('Setup') and \
               self.runParserObj.runparameters.data['RunParameters']['Setup'].has_key('RunMode'):
                return self.runParserObj.runparameters.data['RunParameters']['Setup']['RunMode']
            else:
                raise RuntimeError("not able to guess run mode from RunParameters.xml, parsing problem or new version of software are the likely causes")
        else:
            raise RuntimeError("runParseObj not available")



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
                                




    def demultiplex_run(self):
        """
        Demultiplex a HiSeq run:
            - find the samplesheet
            - make a local copy of the samplesheet and name it SampleSheet.csv
            - create multiple SampleSheets in case at least one lane have multiple indexes lengths
            - run bcl2fastq conversion
        """


        ssname   = self._get_samplesheet()
        ssparser = SampleSheetParser(ssname)
        #Copy the original samplesheet locally. Copy again if already done as there might have been changes to the samplesheet
        try:
            shutil.copy(ssname, self.run_dir)
            ssname = os.path.join(self.run_dir, os.path.split(ssname)[1])
        except:
            raise RuntimeError("unable to copy file {} to destination {}".format(ssname, self.run_dir))
        
        #this sample sheet has been created by the LIMS and copied by a sequencing operator. It is not ready
        #to be used it needs some editing
        #this will contain the samplesheet with all the renaiming to be sued with bcl2fastq-2.17
        samplesheet_dest = os.path.join(self.run_dir, "SampleSheet.csv")
        #check that the samplesheet is not already present. In this case go the next step
        if not os.path.exists(samplesheet_dest):
            try:
                with open(samplesheet_dest, 'wb') as fcd:
                    fcd.write(_generate_clean_samplesheet(ssparser))
            except Exception as e:
                logger.error(e.text)
                return False
            logger.info(("Created SampleSheet.csv for Flowcell {} in {} ".format(self.id, samplesheet_dest)))
        ##SampleSheet.csv generated
        ##when demultiplexing SampleSheet.csv is the one I need to use
        self.runParserObj.samplesheet  = SampleSheetParser(os.path.join(self.run_dir, "SampleSheet.csv"))
        #now geenrate the base masks per lane and decide how to demultiplex
        per_lane_base_masks = self._generate_per_lane_base_mask()
        max_different_base_masks =  max([len(per_lane_base_masks[base_masks]) for base_masks in per_lane_base_masks])
        #if max_different is one, then I have a simple config and I can run a single command. Otherwirse I need to run multiples instances
        #extract lanes with a single base masks
        simple_lanes  = {}
        complex_lanes = {}
        for lane in per_lane_base_masks:
            if len(per_lane_base_masks[lane]) == 1:
                simple_lanes[lane] = per_lane_base_masks[lane]
            else:
                complex_lanes[lane] = per_lane_base_masks[lane]
        #simple lanes contains the lanes such that there is more than one base mask
        bcl2fastq_commands = []
        bcl2fastq_command_num = 0
        if len(simple_lanes) > 0:
            bcl2fastq_commands.append(self._generate_bcl2fastq_command(simple_lanes, True, bcl2fastq_command_num))
            bcl2fastq_command_num += 1
        #compute the different masks, there will be one bcl2fastq command per mask
        base_masks_complex = [complex_lanes[base_masks].keys() for base_masks in complex_lanes]
        different_masks    = list(set([item for sublist in base_masks_complex for item in sublist]))
        for mask in different_masks:
            base_masks_complex_to_demux = {}
            for lane in complex_lanes:
                if complex_lanes[lane].has_key(mask):
                    base_masks_complex_to_demux[lane] = {}
                    base_masks_complex_to_demux[lane][mask] = complex_lanes[lane][mask]
            #at this point base_masks_complex_to_demux contains only a base mask for lane. I can build the command
            bcl2fastq_commands.append(self._generate_bcl2fastq_command(base_masks_complex_to_demux, True, bcl2fastq_command_num))
            bcl2fastq_command_num += 1
        #now bcl2fastq_commands contains all command to be executed. They can be executed in parallel, however run only one per time in order to avoid to overload the machine
        with chdir(self.run_dir):
            # create Demultiplexing dir, in this way the status of this run will became IN_PROGRESS
            if not os.path.exists("Demultiplexing"):
                os.makedirs("Demultiplexing")
            execution = 0
            for bcl2fastq_command in bcl2fastq_commands:
                misc.call_external_command(bcl2fastq_command, with_log_files=True, prefix="demux_{}".format(execution))
                execution += 1



    def _generate_bcl2fastq_command(self, base_masks, strict=True, suffix=0):
        """
        Generates the command to demultiplex with the given base_masks. 
        if strict is set to true demultiplex only lanes in base_masks
        """
        logger.info('Building bcl2fastq command')
        cl = [self.CONFIG.get('bcl2fastq')['bin']]
        if self.CONFIG.get('bcl2fastq').has_key('options'):
            cl_options = self.CONFIG['bcl2fastq']['options']
            # Append all options that appear in the configuration file to the main command.
            for option in cl_options:
                if isinstance(option, dict):
                    opt, val = option.items()[0]
                    #skip output-dir has I might need more than one
                    if "output-dir" not in opt:
                        cl.extend(['--{}'.format(opt), str(val)])
                else:
                    cl.append('--{}'.format(option))
        #now add the base_mask for each lane
        tiles = []
        samplesheetMaskSpecific = os.path.join(os.path.join(self.run_dir, "SampleSheet_{}.csv".format(suffix)))
        output_dir = "Demultiplexing_{}".format(suffix)
        cl.extend(["--output-dir", output_dir])
        
        with open(samplesheetMaskSpecific, 'wb') as ssms:
            ssms.write("[Header]\n")
            ssms.write("[Data]\n")
            ssms.write(",".join(self.runParserObj.samplesheet.datafields))
            ssms.write("\n")
            for lane in sorted(base_masks):
                #iterate thorugh each lane and add the correct --use-bases-mask for that lane
                #there is a single basemask for each lane, I checked it a couple of lines above
                base_mask = [base_masks[lane][bm]['base_mask'] for bm in base_masks[lane]][0] # get the base_mask
                base_mask_expr = "{}:".format(lane) + ",".join(base_mask)
                cl.extend(["--use-bases-mask", base_mask_expr])
                if strict:
                    tiles.extend(["s_{}".format(lane)])
                #these are all the samples that need to be demux with this samplemask in this lane
                samples   = [base_masks[lane][bm]['data'] for bm in base_masks[lane]][0]
                for sample in samples:
                    for field in self.runParserObj.samplesheet.datafields:
                        if field == "index" and "NoIndex" in sample[field]:
                            ssms.write(",") # this is emtpy due to NoIndex issue
                        else:
                            ssms.write("{},".format(sample[field]))
                    ssms.write("\n")
            if strict:
                cl.extend(["--tiles", ",".join(tiles) ])
        cl.extend(["--sample-sheet", samplesheetMaskSpecific])
        logger.info(("BCL to FASTQ command built {} ".format(" ".join(cl))))
        return cl


    def check_run_status(self):
        """
        This function checks the status of a run while in progress.
        In the case of HiSeq check that all demux have been done and in that case perform aggregation
        """
        run_dir    =  self.run_dir
        dex_status =  self.get_run_status()
        #in this case I have already finished all demux jobs and I have aggregate all stasts unded Demultiplexing
        if  dex_status == 'COMPLETED':
            return
        #otherwise check the status of running demux
        #collect all samplesheets generated before
        samplesheets =  glob.glob(os.path.join(run_dir, "*_*.csv"))
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


    def _aggregate_demux_results(self):
        """
        This function aggregates the results from different demultiplexing steps
        """
        run_dir      =  self.run_dir
        demux_folder =  os.path.join(self.run_dir , self.demux_dir)
        samplesheets =  glob.glob(os.path.join(run_dir, "*_*.csv"))
        
        per_lane_base_masks = self._generate_per_lane_base_mask()
        max_different_base_masks =  max([len(per_lane_base_masks[base_masks]) for base_masks in per_lane_base_masks])
        simple_lanes  = {}
        complex_lanes = {}
        for lane in per_lane_base_masks:
            if len(per_lane_base_masks[lane]) == 1:
                simple_lanes[lane] = per_lane_base_masks[lane]
            else:
                complex_lanes[lane] = per_lane_base_masks[lane]
        #complex lanes contains the lanes such that there is more than one base mask
        #for simple lanes undetermined stats will be copied
        if len(complex_lanes) == 0:
            #it means that each lane had only one type of index size, so no need to do super tricky stuff
            print "implement the simple case"
            return
        
        
        html_reports_lane        = []
        html_reports_laneBarcode = []
        for samplesheet in samplesheets:
            demux_id = os.path.splitext(os.path.split(samplesheet)[1])[0].split("_")[1]
            #demux folder is
            demux_id_folder  = os.path.join(run_dir, "Demultiplexing_{}".format(demux_id))
            html_report_lane = os.path.join(run_dir, "Demultiplexing_{}".format(demux_id), "Reports", "html",self.flowcell_id, "all", "all", "all", "lane.html")
            if os.path.exists(html_report_lane):
                html_reports_lane.append(html_report_lane)
            else:
                raise RuntimeError("Not able to find html report {}: possible cause is problem in demultiplexing".format(html_report_lane))
            
            html_report_laneBarcode = os.path.join(run_dir, "Demultiplexing_{}".format(demux_id), "Reports", "html",self.flowcell_id, "all", "all", "all", "laneBarcode.html")
            if os.path.exists(html_report_laneBarcode):
                html_reports_laneBarcode.append(html_report_laneBarcode)
            else:
                raise RuntimeError("Not able to find html report {}: possible cause is problem in demultiplexing".format(html_report_laneBarcode))
            #aggregate fastq
            projects = [project for project in  os.listdir(demux_id_folder) if os.path.isdir(os.path.join(demux_id_folder,project))]
            for project in projects:
                if project in "Reports" or project in "Stats":
                    continue
                project_source = os.path.join(demux_id_folder, project)
                project_dest   = os.path.join(demux_folder, project)
                if not os.path.exists(project_dest):
                    os.makedirs(project_dest)
                #now parse all samples in this project NB: there might be projects that have been demux
                #with different index lenghts, but NO the same sample
                samples = [sample for sample in  os.listdir(project_source) if os.path.isdir(os.path.join(project_source,sample))]
                for sample in samples:
                    sample_source = os.path.join(project_source,sample)
                    sample_dest   = os.path.join(project_dest,sample)
                    if os.path.exists(sample_dest):
                        raise RuntimeError("Sample {} of project {} found twice while aggregating results".format(sample, project))
                    os.makedirs(sample_dest)
                    #now soflink the fastq.gz
                    fastqfiles =  glob.glob(os.path.join(sample_source, "*.fastq*"))
                    for fastqfile in fastqfiles:
                        os.symlink(fastqfile, os.path.join(sample_dest,os.path.split(fastqfile)[1]))
        #now copy fastq files for undetermined (for simple lanes only)
        for lane in simple_lanes.keys():
            undetermined_fastq_files = glob.glob(os.path.join(run_dir, "Demultiplex_0", "Undetermined_S0_*.fastq*")) #should contain only simple lanes undetermined
            for fastqfile in undetermined_fastq_files:
                os.symlink(fastqfile, os.path.join(demux_folder,os.path.split(fastqfile)[1]))
        #now create the html reports
        #start with the lane
        html_report_lane_parser = None
        for next_html_report_lane in html_reports_lane:
            if html_report_lane_parser is None:
                html_report_lane_parser = LaneBarcodeParser(next_html_report_lane)
            else:
                lanesInReport = [Lane['Lane'] for Lane in html_report_lane_parser.sample_data]
                next_html_report_lane_parser = LaneBarcodeParser(next_html_report_lane)
                for entry in next_html_report_lane_parser.sample_data:
                    if not entry["Lane"] in lanesInReport:
                        #if this is a new lane not included before
                        html_report_lane_parser.sample_data.push(entry)
        # now all lanes have been inserted
        for entry in html_report_lane_parser.sample_data:
            if entry['Lane'] in complex_lanes.keys():
                entry['% Perfectbarcode']      = None
                entry['% One mismatchbarcode'] = None
        #now add lanes not present in this demux
        #now I can create the new lane.html
        new_html_report_lane_dir = _create_folder_structure(demux_folder, ["Reports", "html", self.flowcell_id, "all", "all", "all"])
        new_html_report_lane = os.path.join(new_html_report_lane_dir, "lane.html")
        _generate_lane_html(new_html_report_lane, html_report_lane_parser)

        #now generate the laneBarcode
        html_report_laneBarcode_parser = None
        for next_html_report_laneBarcode in html_reports_laneBarcode:
            if html_report_laneBarcode_parser is None:
                html_report_laneBarcode_parser = LaneBarcodeParser(next_html_report_laneBarcode)
            else:
                #no need to check samples occuring in more than one file has I would have spot it while softlinking
                next_html_report_laneBarcode_parser = LaneBarcodeParser(next_html_report_laneBarcode)
                for entry in next_html_report_laneBarcode_parser.sample_data:
                    html_report_laneBarcode_parser.sample_data.append(entry)

        positions_to_delete = [] #find all position that contain default as poriject nameand do not belong to a simple lane
        current_pos = 0
        for entry in html_report_laneBarcode_parser.sample_data:
            if  entry['Lane'] in complex_lanes.keys() and entry['Project'] in "default":
                positions_to_delete = [current_pos] +  positions_to_delete # build the array in this way so that I can delete the elements without messing with the offsets
            current_pos += 1
        for position in positions_to_delete:
            del html_report_laneBarcode_parser.sample_data[position]
        #now generate the new report for laneBarcode.html
        new_html_report_laneBarcode = os.path.join(new_html_report_lane_dir, "laneBarcode.html")
        _generate_lane_html(new_html_report_laneBarcode, html_report_laneBarcode_parser)

        #now create the DemultiplexingStats.xml (empty it is here only to say thay demux is done)
        DemultiplexingStats_xml_dir = _create_folder_structure(demux_folder, ["Stats"])
        #now the run is formally COMPLETED
        open(os.path.join(DemultiplexingStats_xml_dir, "DemultiplexingStats.xml"), 'a').close()



def _create_folder_structure(root, dirs):
    """
    creates a fodler stucture rooted in root usinf all dirs listed in dirs (a list)
    returns the path to the deepest directory
    """
    path=root
    for dir in dirs:
        path = os.path.join(path, dir)
        if not os.path.exists(path):
            os.makedirs(path)
    return path


def _generate_clean_samplesheet(ssparser):
    """
    Will generate a 'clean' samplesheet, for bcl2fastq2.17
    """
    output=""
    #Header
    output+="[Header]{}".format(os.linesep)
    for field in ssparser.header:
        output+="{},{}".format(field.rstrip(), ssparser.header[field].rstrip())
        output+=os.linesep
    
    
    #now parse the data section
    data = []

    for line in ssparser.data:
        entry = {}
        for field, value in line.iteritems():
            if 'SampleID' in field :
                entry[_data_filed_conversion(field)] ='Sample_{}'.format(value)
                entry['Sample_Name'] = value
            elif "Index" in field:
                #in this case we need to distinguish between single and dual index
                entry[_data_filed_conversion(field)] = value.split("-")[0]
                if len(value.split("-")) == 2:
                    entry['index2'] = value.split("-")[1]
                else:
                    entry['index2'] = ""
            else:
                entry[_data_filed_conversion(field)] = value
        data.append(entry)

    fields_to_output = ['Lane', 'Sample_ID', 'Sample_Name', 'index', 'index2', 'Sample_Project']
    #now create the new SampleSheet data section
    output+="[Data]{}".format(os.linesep)
    for field in ssparser.datafields:
        new_field = _data_filed_conversion(field)
        if new_field not in fields_to_output:
            fields_to_output.append(new_field)
    output+=",".join(fields_to_output)
    output+=os.linesep
    #now process each data entry and output it
    for entry in data:
        line = []
        for field in fields_to_output:
            line.append(entry[field])
        output+=",".join(line)
        output+=os.linesep

    return output



def _data_filed_conversion(field):
    """
    converts fields in the sample sheet generated by the LIMS in fields that can be used by bcl2fastq2.17
    """
    datafieldsConversion = {'FCID': 'FCID',
                            'Lane': 'Lane',
                           'SampleID' : 'Sample_ID',
                           'SampleRef': 'SampleRef',
                           'Index' : 'index',
                           'Description': 'Description',
                           'Control': 'Control',
                           'Recipe': 'Recipe',
                           'Operator': 'Operator',
                           'SampleProject' : 'Sample_Project'
                           }
    if field in datafieldsConversion:
        return datafieldsConversion[field]
    else:
        raise RuntimeError("field {} not expected in SampleSheet".format(field))



def _generate_lane_html(html_file, html_report_lane_parser):
    with open(html_file, "w") as html:
        #HEADER
        html.write("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n")
        html.write("<html xmlns:bcl2fastq>\n")
        html.write("<link rel=\"stylesheet\" href=\"../../../../Report.css\" type=\"text/css\">\n")
        html.write("<body>\n")
        html.write("<table width=\"100%\"><tr>\n")
        html.write("<td><p><p>C6L1WANXX /\n")
        html.write("        [all projects] /\n")
        html.write("        [all samples] /\n")
        html.write("        [all barcodes]</p></p></td>\n")
        html.write("<td><p align=\"right\"><a href=\"../../../../FAKE/all/all/all/laneBarcode.html\">show barcodes</a></p></td>\n")
        html.write("</tr></table>\n")
        #FLOWCELL SUMMARY TABLE
        html.write("<h2>Flowcell Summary</h2>\n")
        html.write("<table border=\"1\" ID=\"ReportTable\">\n")
        html.write("<tr>\n")
        keys = html_report_lane_parser.flowcell_data.keys()
        for key in keys:
            html.write("<th>{}</th>\n".format(key))
        html.write("</tr>\n")
        html.write("<tr>\n")
        for key in keys:
            html.write("<td>{}</td>\n".format(html_report_lane_parser.flowcell_data[key]))
        html.write("</tr>\n")
        html.write("</table>\n")
        #LANE SUMMARY TABLE
        html.write("<h2>Lane Summary</h2>\n")
        html.write("<table border=\"1\" ID=\"ReportTable\">\n")
        html.write("<tr>\n")
        keys = html_report_lane_parser.sample_data[0].keys()
        for key in keys:
            html.write("<th>{}</th>\n".format(key))
        html.write("</tr>\n")
        
        for sample in html_report_lane_parser.sample_data:
            html.write("<tr>\n")
            for key in keys:
                html.write("<td>{}</td>\n".format(sample[key]))
            html.write("</tr>\n")
        html.write("</table>\n")
        #FOOTER
        html.write("<p></p>\n")
        html.write("</body>\n")
        html.write("</html>\n")









