from taca.illumina.Standard_Runs import Standard_Run


class NextSeq_Run(Standard_Run):
    def __init__(self, run_dir, software, configuration):
        super(Standard_Run, self).__init__( run_dir, software, configuration)
        self._set_sequencer_type()
        self._set_run_type()
        # NextSeq2000 has a different FC ID pattern that ID contains the first letter for position
        if "VH" in self.instrument:
            self.flowcell_id = self.position + self.flowcell_id
        self._copy_samplesheet()

    def _set_sequencer_type(self):
        self.sequencer_type = "NextSeq"

    def _set_run_type(self):
        self.run_type = "NGI-RUN"
