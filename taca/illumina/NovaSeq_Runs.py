from taca.illumina.HiSeqX_Runs import HiSeqX_Run

import logging

logger = logging.getLogger(__name__)


class NovaSeq_Run(HiSeqX_Run):
    def __init__(self,  run_dir, samplesheet_folders):
        super(NovaSeq_Run, self).__init__( run_dir, samplesheet_folders)
        self._set_sequencer_type()
        self._set_run_type()

    def _set_sequencer_type(self):
        self.sequencer_type = "NovaSeq"

    def _set_run_type(self):
        self.run_type = "NGI-RUN"
