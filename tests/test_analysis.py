#!/usr/bin/env python
import json
import os
import shutil
import tempfile
import unittest
from unittest import mock

from taca.analysis import analysis as an
from taca.utils import config

CONFIG = config.load_yaml_config("data/taca_test_cfg.yaml")


class TestAnalysis(unittest.TestCase):
    """Tests for the Analysis functions."""

    @classmethod
    def setUpClass(self):
        """Creates the following directory tree for testing purposes:

        tmp/
        |__ 141124_ST-COMPLETED_01_AFCIDXX
        |   |__ RunInfo.xml
        |   |__ Demultiplexing
        |   |   |__ Undetermined_S0_L001_R1_001.fastq.gz
        |   |   |__ Stats
        |   |       |__ DemultiplexingStats.xml
        |   |__ RTAComplete.txt
        |   |__ SampleSheet.csv
        """
        self.tmp_dir = os.path.join(tempfile.mkdtemp(), "tmp")
        self.completed = os.path.join(self.tmp_dir, "141124_ST-COMPLETED1_01_AFCIDXX")

        # Create runs directory structure
        os.makedirs(self.tmp_dir)
        os.makedirs(os.path.join(self.completed, "Demultiplexing", "Stats"))

        # Set up files
        open(os.path.join(self.completed, "RTAComplete.txt"), "w").close()
        shutil.copy(
            "data/samplesheet.csv", os.path.join(self.completed, "SampleSheet.csv")
        )
        open(
            os.path.join(
                self.completed, "Demultiplexing", "Stats", "DemultiplexingStats.xml"
            ),
            "w",
        ).close()
        open(
            os.path.join(
                self.completed, "Demultiplexing", "Undetermined_S0_L001_R1_001.fastq.gz"
            ),
            "w",
        ).close()
        with open(
            os.path.join(self.completed, "Demultiplexing", "Stats", "Stats.json"), "w"
        ) as stats_json:
            json.dump({"silly": 1}, stats_json)
        shutil.copy("data/RunInfo.xml", self.completed)
        shutil.copy("data/runParameters.xml", self.completed)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.tmp_dir)

    def test_get_runObj_miseq(self):
        """Return MiSeq run object."""
        miseq_run = os.path.join(self.tmp_dir, "141124_ST-MISEQ1_01_AFCIDXX")
        os.mkdir(miseq_run)
        shutil.copy(
            "data/runParameters_miseq.xml", os.path.join(miseq_run, "runParameters.xml")
        )
        got_miseq_run = an.get_runObj(miseq_run)
        self.assertEqual(got_miseq_run.sequencer_type, "MiSeq")

    def test_get_runObj_nextseq(self):
        """Return NextSeq run object."""
        nextseq_run = os.path.join(self.tmp_dir, "141124_ST-NEXTSEQ1_01_AFCIDXX")
        os.mkdir(nextseq_run)
        shutil.copy(
            "data/runParameters_nextseq.xml",
            os.path.join(nextseq_run, "runParameters.xml"),
        )
        got_nextseq_run = an.get_runObj(nextseq_run)
        self.assertEqual(got_nextseq_run.sequencer_type, "NextSeq")

    def test_get_runObj_novaseq(self):
        """Return NovaSeq run object."""
        novaseq_run = os.path.join(self.tmp_dir, "141124_ST-NOVASEQ1_01_AFCIDXX")
        os.mkdir(novaseq_run)
        shutil.copy(
            "data/runParameters_novaseq.xml",
            os.path.join(novaseq_run, "RunParameters.xml"),
        )
        got_novaseq_run = an.get_runObj(novaseq_run)
        self.assertEqual(got_novaseq_run.sequencer_type, "NovaSeq")

    @mock.patch("taca.analysis.analysis.get_runObj")
    @mock.patch("taca.analysis.analysis._upload_to_statusdb")
    def test_upload_to_statusdb(self, mock_upload_to_statusdb, mock_get_runobj):
        """Get run object and initiate upload to statusdb."""
        mock_get_runobj.return_value = "Standard_run_object"
        an.upload_to_statusdb(self.completed)
        mock_upload_to_statusdb.assert_called_once_with("Stan_run_object")

    @mock.patch("taca.analysis.analysis.statusdb")
    def test__upload_to_statusdb(self, mock_statusdb):
        """Upload to statusdb."""
        run = os.path.join(self.tmp_dir, "141124_ST-NOINDEX1_01_AFCIDYX")
        os.mkdir(run)
        shutil.copy(
            "data/runParameters_minimal.xml", os.path.join(run, "runParameters.xml")
        )
        demux_dir = os.path.join(run, "Demultiplexing", "Stats")
        os.makedirs(demux_dir)
        shutil.copy("data/DemuxSummaryF1L1.txt", demux_dir)
        reports_dir = os.path.join(
            run, "Demultiplexing", "Reports", "html", "FCIDYX", "all", "all", "all"
        )
        os.makedirs(reports_dir)
        shutil.copy("data/laneBarcode.html", (reports_dir))
        shutil.copy("data/lane.html", (reports_dir))
        noindex_run = an.get_runObj(run)
        an._upload_to_statusdb(noindex_run)
        mock_statusdb.update_doc.assert_called_once()

    @mock.patch("taca.analysis.analysis.Standard_Run.transfer_run")
    def test_transfer_run(self, mock_transfer_run):
        """Transfer run to Uppmax."""
        run_dir = self.completed
        an.transfer_run(run_dir)
        mock_transfer_run.assert_called_once_with(
            "nosync/data/transfer.tsv", "some_user@some_email.com"
        )

    @mock.patch("taca.analysis.analysis.RsyncAgent.transfer")
    @mock.patch("taca.analysis.analysis.subprocess.call")
    @mock.patch("taca.analysis.analysis.os.remove")
    @mock.patch("taca.analysis.analysis.open")
    def test_transfer_runfolder(
        self, mock_open, mock_remove, mock_subprocess_call, mock_transfer
    ):
        """Transfer runfolder to uppmax."""
        run_dir = self.completed
        pid = "P1775"
        exclude_lane = ""
        an.transfer_runfolder(run_dir, pid, exclude_lane)
        mock_subprocess_call.assert_called()
        mock_transfer.assert_called()

    def test_extract_project_samplesheet(self):
        """Extract project specific lines from sample sheet."""
        sample_sheet = "data/samplesheet.csv"
        pid = "P1775"
        samplesheet_content = an.extract_project_samplesheet(sample_sheet, pid)
        expected_samplesheet_content = """Lane,SampleID,SampleName,SamplePlate,SampleWell,index,Project
1,Sample_P1775_147,P1775_147,FCB_150423,1:1,GAATTCGT,J_Lundeberg_14_24
"""
        self.assertEqual(samplesheet_content, expected_samplesheet_content)

    @mock.patch("taca.analysis.analysis.NovaSeq_Run.get_run_status")
    @mock.patch("taca.analysis.analysis._upload_to_statusdb")
    def test_run_preprocessing_sequencing(
        self, mock_upload_to_statusdb, mock_get_run_status
    ):
        """Run preprocess run still sequencing."""
        run = self.completed
        mock_get_run_status.return_value = "SEQUENCING"
        an.run_preprocessing(run, force_trasfer=True, statusdb=True)
        mock_upload_to_statusdb.assert_called_once()

    @mock.patch("taca.analysis.analysis.NovaSeq_Run.get_run_status")
    @mock.patch("taca.analysis.analysis._upload_to_statusdb")
    @mock.patch("taca.analysis.analysis.NovaSeq_Run.demultiplex_run")
    def test_run_preprocessing_to_start(
        self, mock_demultiplex_run, mock_upload_to_statusdb, mock_get_run_status
    ):
        """Run preprocessing start demux."""
        run = self.completed
        mock_get_run_status.return_value = "TO_START"
        an.run_preprocessing(run, force_trasfer=True, statusdb=True)
        mock_upload_to_statusdb.assert_called_once()
        mock_demultiplex_run.assert_called_once()

    @mock.patch("taca.analysis.analysis.NovaSeq_Run.get_run_status")
    @mock.patch("taca.analysis.analysis._upload_to_statusdb")
    @mock.patch("taca.analysis.analysis.NovaSeq_Run.check_run_status")
    def test_run_preprocessing_in_progress(
        self, mock_check_run_status, mock_upload_to_statusdb, mock_get_run_status
    ):
        """Run preprocessing demux in progress."""
        run = self.completed
        mock_get_run_status.return_value = "IN_PROGRESS"
        an.run_preprocessing(run, force_trasfer=True, statusdb=True)
        mock_upload_to_statusdb.assert_called_once()
        mock_check_run_status.assert_called_once()

    @mock.patch("taca.analysis.analysis.NovaSeq_Run.get_run_status")
    @mock.patch("taca.analysis.analysis._upload_to_statusdb")
    @mock.patch("taca.analysis.analysis.NovaSeq_Run.send_mail")
    @mock.patch("taca.analysis.analysis.NovaSeq_Run.transfer_run")
    @mock.patch("taca.analysis.analysis.os.mkdir")
    @mock.patch("taca.analysis.analysis.copyfile")
    def test_run_preprocessing_completed(
        self,
        mock_copy,
        mock_mkdir,
        mock_transfer_run,
        mock_send_mail,
        mock_upload_to_statusdb,
        mock_get_run_status,
    ):
        """Run preprocessing demux completed."""
        run = self.completed
        mock_get_run_status.return_value = "COMPLETED"
        an.run_preprocessing(run, force_trasfer=True, statusdb=True)
        mock_upload_to_statusdb.assert_called_once()
        message = "The run 141124_ST-COMPLETED1_01_AFCIDXX has been demultiplexed.\n                The Run will be transferred to the analysis cluster for further analysis.\n\n             \
   The run is available at : https://genomics-status.scilifelab.se/flowcells/141124_ST-COMPLETED1_01_AFCIDXX\n\n                "
        mock_send_mail.assert_called_once_with(message, rcp="some_user@some_email.com")
        mock_transfer_run.assert_called_once_with(
            "data/transfer.tsv", "some_user@some_email.com"
        )
