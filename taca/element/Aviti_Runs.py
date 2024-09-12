from taca.element.Element_Runs import Run
from taca.utils.config import CONFIG
from taca.utils.statusdb import ElementRunsConnection


class Aviti_Run(Run):
    def __init__(self, run_dir, configuration):
        super().__init__(run_dir, configuration)
        self.sequencer_type = "Aviti"
        self.db = ElementRunsConnection(CONFIG["statusdb"], dbname="element_runs")

    def update_statusdb(self):
        doc_obj = self.to_doc_obj()
        self.db.upload_to_statusdb(doc_obj)

    def construct_NGI_run_id(self):
        pass

    def to_doc_obj(self):
        doc_obj = {
            "name": self.run_name,
            "run_status": self.run_status,
            "run_id": self.run_id,
            "run_dir": self.run_dir,
            "run_type": self.run_type,
            "sequencer_type": self.sequencer_type,
            "samples": self.samples,
            "demux": self.demux,
            "date": self.date,
            "flowcell": self.flowcell,
        }
        return doc_obj
