from taca.illumina.HiSeqX_Runs import HiSeqX_Run


class NextSeq_Run(HiSeqX_Run):

    def __init__(self,  run_dir, samplesheet_folders):
        super(HiSeqX_Run, self).__init__(run_dir, samplesheet_folders)
        self._set_sequencer_type()
        self._set_run_type()
        self.flowcell_id = self.position + self.flowcell_id
        self._copy_samplesheet()

    def _set_sequencer_type(self):
        self.sequencer_type = "NextSeq"

    def _set_run_type(self):
        self.run_type = "NGI-RUN"
