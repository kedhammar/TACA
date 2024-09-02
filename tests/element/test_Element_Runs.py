import os
import tempfile

import pytest

from taca.element import Element_Runs as to_test


def create_aviti_run_dir(
    tmp: tempfile.TemporaryDirectory,
    run_name: str = "20240716_AV242106_testrun",
    nosync: bool = False,
    run_finished: bool = True,
    sync_finished: bool = True,
) -> str:
    # Create run dir
    if nosync:
        run_path = f"{tmp.name}/ngi_data/sequencing/AV242106/nosync/{run_name}"
    else:
        run_path = f"{tmp.name}/ngi_data/sequencing/AV242106/{run_name}"
    os.mkdir(run_path)

    # Create files
    if run_finished:
        open(f"{run_path}/AvitiRunStats.json", "w").close()
        open(f"{run_path}/RunManifest.csv", "w").close()
        open(f"{run_path}/RunManifest.json", "w").close()
        open(f"{run_path}/RunParameters.json", "w").close()
        open(f"{run_path}/RunUploaded.json", "w").close()

    if sync_finished:
        open(f"{run_path}/.sync_finished", "w").close()

    return run_path


class TestRun:
    @pytest.fixture(autouse=True)
    def setup(self, create_dirs: pytest.fixture):
        self.tmp: tempfile.TemporaryDirectory = create_dirs
        self.run_path = create_aviti_run_dir(self.tmp)
        self.run = to_test.Run(self.run_path, {})

    def test_init(self):
        assert self.run.run_dir == self.run_path

    def test_check_sequencing_status(self):
        assert False

    def test_get_demultiplexing_status(self):
        assert False

    def test_manifest_exists(self):
        assert False

    def test_generate_demux_command(self):
        assert False

    def test_start_demux(self):
        assert False

    def test_is_transferred(self):
        assert False

    def test_parse_rundir(self):
        assert False
