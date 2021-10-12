from taca.nanopore.nanopore import Nanopore

class PromethION(Nanopore):
    """PromethION run"""
    def __init__(self, run_dir):
        super(PromethION, self).__init__(run_dir)
    