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
from taca.illumina.HiSeq_Runs import HiSeq_Run
from taca.utils import misc
from flowcell_parser.classes import RunParametersParser, SampleSheetParser, RunParser, LaneBarcodeParser, DemuxSummaryParser


import logging

logger = logging.getLogger(__name__)

class MiSeq_Run(HiSeq_Run):

    def __init__(self,  path_to_run, configuration):
        super(MiSeq_Run, self).__init__( path_to_run, configuration)
        self._set_sequencer_type()


    def _set_sequencer_type(self):
        self.sequencer_type = "MiSeq"
    
    
    def _sequencer_type(self):
            return "MiSeq"




    def _get_samplesheet(self):
        """
        Locate and parse the samplesheet for a run.
        In MiSeq case this is located in FC_DIR/Data/Intensities/BaseCalls/SampleSheet.csv
        """
        ssname = os.path.join(self.run_dir, 'Data', 'Intensities', 'BaseCalls','SampleSheet.csv')
        if os.path.exists(ssname):
            #if exists parse the SampleSheet
            ssparser = SampleSheetParser(ssname)
            if ssparser.header['Description'] == "Production" or ssparser.header['Description'] == "Application":
                # if Description is either Production or Application we need to process it
                return ssname
            else:
                #otherwise foget about it
                self.archive_run(os.path.join(os.path.split(self.run_dir)[0], "nosync"))
                return None
        else:
            #In this case move the run to NoSync and forget about it
            self.archive_run(os.path.join(os.path.split(self.run_dir)[0], "nosync"))
            return None

                                


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
                if 'Sample_ID' in field:
                    entry[field] ='Sample_{}'.format(value)
                elif 'Sample_Project' in field:
                    entry[field] = value.replace(".", "_")
                else:
                    entry[field] = value
            if 'Lane' not in entry:
                entry['Lane'] = '1'
            data.append(entry)

        fields_to_output = ['Lane', 'Sample_ID', 'Sample_Name', 'index', 'index2', 'Sample_Project']
        #now create the new SampleSheet data section
        output+="[Data]{}".format(os.linesep)
        for field in ssparser.datafields:
            if field not in fields_to_output:
                fields_to_output.append(field)
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






