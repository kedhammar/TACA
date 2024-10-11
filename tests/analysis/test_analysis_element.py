from tempfile import TemporaryDirectory
from unittest.mock import patch

from tests.element.test_Element_Runs import create_element_run_dir, get_config


def test_run_preprocessing(create_dirs):
    tmp: TemporaryDirectory = create_dirs

    # Mock config
    config = get_config(tmp)
    mock_config = patch("taca.utils.config.CONFIG", new=config)
    mock_config.start()

    # Mock DB
    mock_db = patch("taca.element.Element_Runs.ElementRunsConnection")
    mock_db.start()

    # Mock subprocess
    mock_subprocess = patch("subprocess.Popen")
    mock_subprocess.start()

    # Create run dir and associated LIMS manifest
    run_dir = create_element_run_dir(tmp=tmp)

    # Import module to test
    from taca.analysis import analysis_element as to_test

    # Test
    to_test.run_preprocessing(run_dir)

    # Stop mocks
    patch.stopall()
