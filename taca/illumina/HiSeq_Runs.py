import os
import re
import csv
import glob
import shutil
import gzip
import operator
import subprocess
from datetime import datetime
from taca.utils.filesystem import chdir, control_fastq_filename
from taca.illumina.Runs import Run
from taca.utils import misc
from flowcell_parser.classes import RunParametersParser, SampleSheetParser, RunParser, LaneBarcodeParser, DemuxSummaryParser


import logging

logger = logging.getLogger(__name__)

class HiSeq_Run(Run):

    def __init__(self,  path_to_run, samplesheet_folders):
        super(HiSeq_Run, self).__init__( path_to_run, samplesheet_folders)
        self._set_sequencer_type()
        self._set_run_type()

    def _set_sequencer_type(self):
        self.sequencer_type = "HiSeq"
    
    def _set_run_type(self):
        self.run_type = "NGI-RUN"

    def _get_run_mode(self): #Old function, not really used but might be usefull in the future
        if self.runParserObj:
            if self.runParserObj.runparameters.data.has_key('RunParameters') and \
               self.runParserObj.runparameters.data['RunParameters'].has_key('Setup') and \
               self.runParserObj.runparameters.data['RunParameters']['Setup'].has_key('RunMode'):
                return self.runParserObj.runparameters.data['RunParameters']['Setup']['RunMode']
            else:
                raise RuntimeError("not able to guess run mode from RunParameters.xml, parsing problem or new version of software are the likely causes")
        else:
            raise RuntimeError("runParseObj not available")

    def demultiplex_run(self):
        """
        Demultiplex a HiSeq run:
            - find the samplesheet
            - make a local copy of the samplesheet and name it SampleSheet.csv
            - create multiple SampleSheets in case at least one lane have multiple indexes lengths
            - run bcl2fastq conversion
        """

        ssname   = self._get_samplesheet()
        if ssname is None:
            return None
        ssparser = SampleSheetParser(ssname)
        #Copy the original samplesheet locally. Copy again if already done as there might have been changes to the samplesheet
        try:
            shutil.copy(ssname, os.path.join(self.run_dir, "{}.csv".format(self.flowcell_id)))
            ssname = os.path.join(self.run_dir, os.path.split(ssname)[1])
        except:
            raise RuntimeError("unable to copy file {} to destination {}".format(ssname, self.run_dir))
        
        #this sample sheet has been created by the LIMS and copied by a sequencing operator. It is not ready
        #to be used it needs some editing
        #this will contain the samplesheet with all the renaiming to be used with bcl2fastq-2.17
        samplesheet_dest = os.path.join(self.run_dir, "SampleSheet.csv")
        #check that the samplesheet is not already present. In this case go the next step
        if os.path.exists(samplesheet_dest):
            logger.info("SampleSheet.csv found ... overwriting it")
        try:
            with open(samplesheet_dest, 'wb') as fcd:
                fcd.write(self._generate_clean_samplesheet(ssparser))
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
                misc.call_external_command_detached(bcl2fastq_command, with_log_files=True, prefix="demux_{}".format(execution))
                execution += 1



    def _generate_bcl2fastq_command(self, base_masks, strict=True, suffix=0, mask_short_adapter_reads=False):
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
                        if field == "index" and "NOINDEX" in sample[field]:
                            ssms.write(",") # this is emtpy due to NoIndex issue
                        else:
                            ssms.write("{},".format(sample[field]))
                    ssms.write("\n")
            if strict:
                cl.extend(["--tiles", ",".join(tiles) ])
        cl.extend(["--sample-sheet", samplesheetMaskSpecific])
        if mask_short_adapter_reads:
            cl.extend(["--mask-short-adapter-reads", "0"])
        
        logger.info(("BCL to FASTQ command built {} ".format(" ".join(cl))))
        return cl



    def _aggregate_demux_results(self):
        """
        This function aggregates the results from different demultiplexing steps
        """
        run_dir      =  self.run_dir
        demux_folder =  os.path.join(self.run_dir , self.demux_dir)
        samplesheets =  glob.glob(os.path.join(run_dir, "*_[0-9].csv")) # a single digit... this hipotesis should hold for a while
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
            demux_folder_tmp_name = "Demultiplexing_0" # in this case this is the only demux dir
            demux_folder_tmp     = os.path.join(run_dir, demux_folder_tmp_name)
            elements = [element for element  in  os.listdir(demux_folder_tmp) ]
            for element in elements:
                if "Stats" not in element: #skip this folder and treat it differently to take into account the NoIndex case
                    source  = os.path.join(demux_folder_tmp, element)
                    dest    = os.path.join(self.run_dir, self.demux_dir, element)
                    os.symlink(source, dest)
            os.makedirs(os.path.join(self.run_dir, "Demultiplexing", "Stats"))
            #now fetch the lanes that have NoIndex
            noIndexLanes = [Sample["Lane"] for Sample in  self.runParserObj.samplesheet.data if "NOINDEX" in Sample["index"]]
            statsFiles = glob.glob(os.path.join(demux_folder_tmp, "Stats", "*" ))
            for source in statsFiles:
                source_name = os.path.split(source)[1]
                if source_name not in ["DemultiplexingStats.xml", "AdapterTrimming.txt", "ConversionStats.xml"]:
                    lane = os.path.splitext(os.path.split(source)[1])[0][-1] #lane
                    if lane not in noIndexLanes:
                        #in this case I can soflink the file here
                        dest    = os.path.join(self.run_dir, self.demux_dir, "Stats", source_name)
                        os.symlink(source, dest)
            #now copy the three last files
            for file in ["DemultiplexingStats.xml", "AdapterTrimming.txt", "ConversionStats.xml"]:
                source = os.path.join(self.run_dir, "Demultiplexing_0", "Stats", file)
                dest   = os.path.join(self.run_dir, "Demultiplexing", "Stats", file)
                os.symlink(source, dest)
            #this is the simple case, Demultiplexing dir is simply a symlink to the only sub-demultiplexing dir
            return True
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
                    #there might be project seqeunced with multiple index lengths
                    os.makedirs(project_dest)
                samples = [sample for sample in  os.listdir(project_source) if os.path.isdir(os.path.join(project_source,sample))]
                for sample in samples:
                    sample_source = os.path.join(project_source,sample)
                    sample_dest   = os.path.join(project_dest,sample)
                    if not os.path.exists(sample_dest):
                        #there should beven be the same sample sequenced with different index length, however a sample might be pooled in several lanes and therefore sequenced using different samplesheets.
                        os.makedirs(sample_dest)
                    #now soflink the fastq.gz
                    fastqfiles =  glob.glob(os.path.join(sample_source, "*.fastq*"))
                    for fastqfile in fastqfiles:
                        os.symlink(fastqfile, os.path.join(sample_dest,os.path.split(fastqfile)[1]))

        #now copy fastq files for undetermined (for simple lanes only)
        for lane in simple_lanes.keys():
            undetermined_fastq_files = glob.glob(os.path.join(run_dir, "Demultiplexing_0", "Undetermined_S0_L00{}*.fastq*".format(lane))) #contains only simple lanes undetermined
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
                        html_report_lane_parser.sample_data.append(entry)
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
        #copy the Undetermined stats for simple lanes
        for lane in simple_lanes.keys():
            DemuxSummaryFiles = glob.glob(os.path.join(run_dir, "Demultiplexing_0", "Stats", "*L{}*txt".format(lane)))
            for DemuxSummaryFile in DemuxSummaryFiles:
                os.symlink(DemuxSummaryFile, os.path.join(demux_folder, "Stats", os.path.split(DemuxSummaryFile)[1]))
        #now the run is formally COMPLETED
        open(os.path.join(DemultiplexingStats_xml_dir, "DemultiplexingStats.xml"), 'a').close()
        return True





    def _generate_clean_samplesheet(self, ssparser):
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
                    entry[_data_filed_conversion(field)] = value.split("-")[0].upper()
                    if len(value.split("-")) == 2:
                        entry['index2'] = value.split("-")[1].upper()
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









