from taca.nanopore.nanopore import Nanopore
from taca.utils.config import CONFIG

class PromethION(Nanopore):
    """PromethION run"""
    def __init__(self, run_dir):
        super(PromethION, self).__init__(run_dir)
        self.transfer_log = CONFIG.get('nanopore_analysis').get('promethion_run').get('transfer').get('transfer_file')
        self.archive_dir = CONFIG.get('nanopore_analysis').get('promethion_run').get('finished_dir')