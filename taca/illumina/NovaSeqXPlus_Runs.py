from taca.illumina.HiSeqX_Runs import HiSeqX_Run


class NovaSeqXPlus_Run(HiSeqX_Run):
    def __init__(self, run_dir, configuration):
        super(NovaSeqXPlus_Run, self).__init__(run_dir, configuration)
        self._set_sequencer_type()
        self._set_run_type()

    def _set_sequencer_type(self):
        self.sequencer_type = "NovaSeqXPlus"

    def _set_run_type(self):
        self.run_type = "NGI-RUN"

    def _current_year(self):
        """Method needed to extract year from rundir name, since year contains 4 digits
        on NovaSeqXPlus while previously it was 2."""
        return self.id[0:4]