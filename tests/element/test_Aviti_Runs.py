import tempfile

import pytest

from taca.element import Aviti_Runs as to_test
from tests.element.test_Element_Runs import create_element_run_dir


class TestAviti_Run:
    def test_init(self, create_dirs: pytest.fixture):
        tmp: tempfile.TemporaryDirectory = create_dirs
        run_dir = create_element_run_dir(tmp)
        run = to_test.Aviti_Run(run_dir, {})
        assert run.run_dir == run_dir
        assert run.sequencer_type == "Aviti"
