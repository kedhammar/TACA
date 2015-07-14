import os
import re
import csv
import glob
import datetime
from hugin.Runs.Runs import Run
from hugin.parser import RunInfoParser
from hugin.parser import  HiSeqSampleSheet
from hugin.parser import RunParametersParser

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
        RP = RunParametersParser()
        return RP.parse(os.path.join(self.path, "runParameters.xml"))['RunMode']





