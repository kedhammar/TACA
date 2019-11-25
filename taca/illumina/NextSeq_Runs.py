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
from taca.illumina.HiSeqX_Runs import HiSeqX_Run
from taca.utils import misc
from flowcell_parser.classes import RunParametersParser, SampleSheetParser, RunParser, LaneBarcodeParser, DemuxSummaryParser


import logging

logger = logging.getLogger(__name__)

class NextSeq_Run(HiSeqX_Run):

    def __init__(self,  path_to_run, configuration):
        # Constructor, it returns a NextSeq object only
        # if the NextSeq run belongs to NGI facility, i.e., contains
        # Application or production in the Description
        super(NextSeq_Run, self).__init__( path_to_run, configuration)
        # In the NextSeq the sample sheet is created by the operator
        # and placed in the run root folder.
        # For now we use the flow cell id to identify the sample sheet
        self.ssname = os.path.join(self.run_dir, self.flowcell_id + ".csv")
        self._set_sequencer_type()
        self._set_run_type()

    def _set_sequencer_type(self):
        self.sequencer_type = "NextSeq"

    def _set_run_type(self):
        if not os.path.exists(self.ssname):
            # Case in which no samplesheet is found, assume it is a non NGI run
            self.run_type = "NON-NGI-RUN"
        else:
            # it SampleSheet exists try to see if it is a NGI-run
            # TODO SampleSheetParser may throw an exception
            ssparser = SampleSheetParser(self.ssname)
            if ssparser.header['Description'] == "Production" or ssparser.header['Description'] == "Applications":
                self.run_type = "NGI-RUN"
            else:
            #otherwise this is a non NGI run
                self.run_type = "NON-NGI-RUN"
