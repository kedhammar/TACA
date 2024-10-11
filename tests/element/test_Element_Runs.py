import os
import tempfile
import zipfile
from unittest import mock

import pytest

from taca.element import Element_Runs as to_test


def get_config(tmp: tempfile.TemporaryDirectory) -> dict:
    config = {
        "element_analysis": {
            "Element": {
                "Aviti": {
                    "manifest_zip_location": f"{tmp.name}/ngi-nas-ns/samplesheets/Aviti",
                    "transfer_log": f"{tmp.name}/log/transfer_aviti.tsv",
                },
            },
        },
        "statusdb": {},
    }
    return config


def create_element_run_dir(
    tmp: tempfile.TemporaryDirectory,
    run_name: str = "20240926_AV242106_A2349523513",
    lims_manifest: bool = True,
    nosync: bool = False,
    run_finished: bool = True,
    sync_finished: bool = True,
    demux_dir: bool = True,
    n_demux_subdirs: int = 1,
    demux_done: bool = True,
    outcome_completed: bool = True,
) -> str:
    """
    Build a run dir for an Element run for test purposes.

    Some file contents are replaced with "MOCK" to shorten them.

        20240926_AV242106_A2349523513
        ├── RunManifest.csv
        ├── RunManifest.json
        ├── RunParameters.json
        ├── RunUploaded.json
        ├── .sync_finished
        ├── Demultiplexing
        ├── Demultiplexing_0
        |   └── RunStats.json
        ├── Demultiplexing_1
        |   └── RunStats.json
        └── ...

    """

    # Create run dir
    if nosync:
        run_path = f"{tmp.name}/ngi_data/sequencing/AV242106/nosync/{run_name}"
    else:
        run_path = f"{tmp.name}/ngi_data/sequencing/AV242106/{run_name}"
    os.mkdir(run_path)

    # Create LIMS manifest
    if lims_manifest:
        manifest_root_name = "AVITI_run_manifest_2349523513_24-1061390_240926_171138_ChristianNatanaelsson"
        manifest_pdir = f"{tmp.name}/ngi-nas-ns/samplesheets/Aviti/2024"

        csv_path = f"{manifest_pdir}/{manifest_root_name}_untrimmed.csv"
        zip_path = f"{manifest_pdir}/{manifest_root_name}.zip"

        with open(csv_path, "w") as stream:
            stream.write("""[RUNVALUES]
KeyName, Value
lims_step_name, "Load to Flowcell (AVITI) v1.0"
lims_step_id, "24-1061390"
manifest_file, "AVITI_run_manifest_2349523513_24-1061390_240926_171138_ChristianNatanaelsson_untrimmed.csv"

[SETTINGS]
SettingName, Value

[SAMPLES]
SampleName,Index1,Index2,Lane,Project,Recipe
P32105_1001,AAAGCATA,,1,I__Adameyko_24_06,50-8-24-49
P32105_1001,CTGCAGCC,,1,I__Adameyko_24_06,50-8-24-49
P32105_1001,GCCTTTAT,,1,I__Adameyko_24_06,50-8-24-49
P32105_1001,TGTAGCGG,,1,I__Adameyko_24_06,50-8-24-49
P32105_1002,ATTGGACG,,1,I__Adameyko_24_06,50-8-24-49
P32105_1002,CAGCTTAC,,1,I__Adameyko_24_06,50-8-24-49
P32105_1002,GGCAAGGA,,1,I__Adameyko_24_06,50-8-24-49
P32105_1002,TCATCCTT,,1,I__Adameyko_24_06,50-8-24-49
P32105_1003,ACGTTACA,,1,I__Adameyko_24_06,50-8-24-49
P32105_1003,CGTAGGTT,,1,I__Adameyko_24_06,50-8-24-49
P32105_1003,GACGACGG,,1,I__Adameyko_24_06,50-8-24-49
P32105_1003,TTACCTAC,,1,I__Adameyko_24_06,50-8-24-49
P32105_1004,ACTTCACT,,1,I__Adameyko_24_06,50-8-24-49
P32105_1004,CGAAGTTG,,1,I__Adameyko_24_06,50-8-24-49
P32105_1004,GAGCACGC,,1,I__Adameyko_24_06,50-8-24-49
P32105_1004,TTCGTGAA,,1,I__Adameyko_24_06,50-8-24-49
PhiX_Adept,ATGTCGCTAG,CTAGCTCGTA,1,Control,0-0
PhiX_Adept,CACAGATCGT,ACGAGAGTCT,1,Control,0-0
PhiX_Adept,GCACATAGTC,GACTACTAGC,1,Control,0-0
PhiX_Adept,TGTGTCGACA,TGTCTGACAG,1,Control,0-0
P32105_1001,AAAGCATA,,2,I__Adameyko_24_06,50-8-24-49
P32105_1001,CTGCAGCC,,2,I__Adameyko_24_06,50-8-24-49
P32105_1001,GCCTTTAT,,2,I__Adameyko_24_06,50-8-24-49
P32105_1001,TGTAGCGG,,2,I__Adameyko_24_06,50-8-24-49
P32105_1002,ATTGGACG,,2,I__Adameyko_24_06,50-8-24-49
P32105_1002,CAGCTTAC,,2,I__Adameyko_24_06,50-8-24-49
P32105_1002,GGCAAGGA,,2,I__Adameyko_24_06,50-8-24-49
P32105_1002,TCATCCTT,,2,I__Adameyko_24_06,50-8-24-49
P32105_1003,ACGTTACA,,2,I__Adameyko_24_06,50-8-24-49
P32105_1003,CGTAGGTT,,2,I__Adameyko_24_06,50-8-24-49
P32105_1003,GACGACGG,,2,I__Adameyko_24_06,50-8-24-49
P32105_1003,TTACCTAC,,2,I__Adameyko_24_06,50-8-24-49
P32105_1004,ACTTCACT,,2,I__Adameyko_24_06,50-8-24-49
P32105_1004,CGAAGTTG,,2,I__Adameyko_24_06,50-8-24-49
P32105_1004,GAGCACGC,,2,I__Adameyko_24_06,50-8-24-49
P32105_1004,TTCGTGAA,,2,I__Adameyko_24_06,50-8-24-49
PhiX_Adept,ATGTCGCTAG,CTAGCTCGTA,2,Control,0-0
PhiX_Adept,CACAGATCGT,ACGAGAGTCT,2,Control,0-0
PhiX_Adept,GCACATAGTC,GACTACTAGC,2,Control,0-0
PhiX_Adept,TGTGTCGACA,TGTCTGACAG,2,Control,0-0
""")

        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
            # Add the CSV file to the zip file
            zipf.write(csv_path, os.path.basename(csv_path))

    # Populate run dir with files and folders
    if run_finished:
        with open(f"{run_path}/RunManifest.json", "w") as stream:
            stream.write("""{
    "KitConfiguration": {
        "MaxCycles": 334,
        "DefaultR1Cycles": 151,
        "DefaultR2Cycles": 151,
        "DefaultI1Cycles": -1,
        "DefaultI2Cycles": -1,
        "MinimumR1Cycles": 5,
        "MinimumR2Cycles": 0,
        "MinimumI1Cycles": 0,
        "MinimumI2Cycles": 0,
        "DefaultI1FastQ": false,
        "DefaultI2FastQ": false,
        "DefaultUMIFastQ": false,
        "DefaultI1Mask": "I1:Y*",
        "DefaultI2Mask": "I2:Y*",
        "DefaultUmiMask": "I1:N*",
        "DefaultR1FastQMask": "R1:Y*N",
        "DefaultR2FastQMask": "R2:Y*N",
        "DefaultI1MaskRead": "I1",
        "DefaultI2MaskRead": "I2",
        "DefaultUmiMaskRead": "I1",
        "DefaultR1FastQMaskRead": "R1",
        "DefaultR2FastQMaskRead": "R2",
        "DefaultR1Adapter": "",
        "DefaultR2Adapter": "",
        "DefaultR1AdapterTrim": false,
        "DefaultR2AdapterTrim": false,
        "DefaultR1AdapterNMask": false,
        "DefaultR2AdapterNMask": false,
        "DefaultR1AdapterMinimumTrimmedLength": 16,
        "DefaultR2AdapterMinimumTrimmedLength": 16,
        "DefaultR1AdapterMinimumStringency": 0.9,
        "DefaultR2AdapterMinimumStringency": 0.9,
        "DefaultR1AdapterMinimumOverlap": 3,
        "DefaultR2AdapterMinimumOverlap": 3,
        "DefaultAdapterTrimType": "Paired-End"
    },
    "RunParameters": {
        "PreparationWorkflow": "Adept",
        "KitConfiguration": "300Cycles",
        "ChemistryVersion": "Cloudbreak",
        "LowDiversity": false,
        "I1Cycles": 8,
        "I2Cycles": 24,
        "R1Cycles": 50,
        "R2Cycles": 49
    },
    "RunValues": {
        "lims_step_id": "24-1061390",
        "lims_step_name": "Load to Flowcell (AVITI) v1.0",
        "manifest_file": "AVITI_run_manifest_2349523513_24-1061390_240926_171138_ChristianNatanaelsson_trimmed.csv"
    },
    "Settings": [
        {
            "Lane": 1,
            "I1MismatchThreshold": 1,
            "I2MismatchThreshold": 1,
            "R1Adapter": [],
            "R2Adapter": [],
            "I1MaskManifest": "I1:N*",
            "I1Mask": [
                {
                    "Read": "I1",
                    "Cycles": []
                }
            ],
            "I1FastQ": false,
            "I2MaskManifest": "I2:N*",
            "I2Mask": [
                {
                    "Read": "I2",
                    "Cycles": []
                }
            ],
            "I2FastQ": false,
            "UmiMaskManifest": "I1:N*",
            "UmiMask": [
                {
                    "Read": "I1",
                    "Cycles": []
                }
            ],
            "UmiFastQ": false,
            "R1FastQMaskManifest": "R1:Y*N",
            "R1FastQMask": [
                {
                    "Read": "R1",
                    "Cycles": "MOCK"
                }
            ],
            "R2FastQMaskManifest": "R2:Y*N",
            "R2FastQMask": [
                {
                    "Read": "R2",
                    "Cycles": "MOCK"
                }
            ],
            "SpikeInAsUnassigned": true,
            "R1AdapterTrim": false,
            "R2AdapterTrim": false,
            "R1AdapterNMask": false,
            "R2AdapterNMask": false,
            "R1AdapterMinimumTrimmedLength": 16,
            "R2AdapterMinimumTrimmedLength": 16,
            "R1AdapterMinimumStringency": 0.9,
            "R2AdapterMinimumStringency": 0.9,
            "R1AdapterMinimumOverlap": 3,
            "R2AdapterMinimumOverlap": 3,
            "AdapterTrimType": "Paired-End"
        },
        {
            "Lane": 2,
            "I1MismatchThreshold": 1,
            "I2MismatchThreshold": 1,
            "R1Adapter": [],
            "R2Adapter": [],
            "I1MaskManifest": "I1:N*",
            "I1Mask": [
                {
                    "Read": "I1",
                    "Cycles": []
                }
            ],
            "I1FastQ": false,
            "I2MaskManifest": "I2:N*",
            "I2Mask": [
                {
                    "Read": "I2",
                    "Cycles": []
                }
            ],
            "I2FastQ": false,
            "UmiMaskManifest": "I1:N*",
            "UmiMask": [
                {
                    "Read": "I1",
                    "Cycles": []
                }
            ],
            "UmiFastQ": false,
            "R1FastQMaskManifest": "R1:Y*N",
            "R1FastQMask": [
                {
                    "Read": "R1",
                    "Cycles": "MOCK"
                }
            ],
            "R2FastQMaskManifest": "R2:Y*N",
            "R2FastQMask": [
                {
                    "Read": "R2",
                    "Cycles": "MOCK"
                }
            ],
            "SpikeInAsUnassigned": true,
            "R1AdapterTrim": false,
            "R2AdapterTrim": false,
            "R1AdapterNMask": false,
            "R2AdapterNMask": false,
            "R1AdapterMinimumTrimmedLength": 16,
            "R2AdapterMinimumTrimmedLength": 16,
            "R1AdapterMinimumStringency": 0.9,
            "R2AdapterMinimumStringency": 0.9,
            "R1AdapterMinimumOverlap": 3,
            "R2AdapterMinimumOverlap": 3,
            "AdapterTrimType": "Paired-End"
        }
    ],
    "Samples": [
        {
            "SampleName": "DefaultSample",
            "SampleNumber": 1,
            "ExternalId": "",
            "Indexes": [
                {
                    "Lane": 1,
                    "Index1": "",
                    "Index2": ""
                },
                {
                    "Lane": 2,
                    "Index1": "",
                    "Index2": ""
                }
            ],
            "CustomMetadata": {},
            "Project": "DefaultProject"
        }
    ]
}
""")
        with open(f"{run_path}/RunParameters.json", "w") as stream:
            stream.write("""{
  "FileVersion": "5.0.0",
  "RunName": "A2349523513",
  "RecipeExecutionID": "rec.9590c80c95fc4eee8b3eb10c31251915",
  "RunID": "seq_66f5837f1ae1a35f10a2e594",
  "RunType": "Sequencing",
  "RunDescription": "",
  "Side": "SideA",
  "FlowcellID": "2349523513",
  "Date": "2024-09-26T16:34:55.978072698Z",
  "InstrumentName": "AV242106",
  "OperatorName": "christian.natanael@scilifelab.se ",
  "RunFolderName": "20240926_AV242106_A2349523513",
  "Tiles": "MOCK",
  "Cycles": {
    "R1": 50,
    "R2": 49,
    "I1": 8,
    "I2": 24
  },
  "ReadOrder": "I1,I2,R1,R2",
  "ThroughputSelection": "High",
  "KitConfiguration": "300Cycles",
  "PreparationWorkflow": "Adept",
  "ChemistryVersion": "Cloudbreak",
  "LowDiversity": false,
  "PlatformVersion": "2.6.2",
  "AnalysisLanes": "1+2",
  "StorageConnectionID": "local:66866355d07c3234c01b67b1",
  "PMGMask": "P1:Y4N*",
  "Consumables": {
    "Flowcell": {
      "SerialNumber": "2349523513",
      "PartNumber": "810-00002",
      "LotNumber": "2405300233",
      "Expiration": "2025-05-31T00:00:00Z",
      "ExpirationStr": "20250531",
      "BarcodeStr": "2349523513,810-00002,2405300233,20250531"
    },
    "SequencingCartridge": {
      "SerialNumber": "24062600390028",
      "PartNumber": "820-00013",
      "LotNumber": "2406260039",
      "Expiration": "2025-05-22T00:00:00Z",
      "ExpirationStr": "20250522",
      "BarcodeStr": "24062600390028,820-00013,2406260039,20250522"
    },
    "Buffer": {
      "SerialNumber": "24062400390041",
      "PartNumber": "820-00002",
      "LotNumber": "2406240039",
      "Expiration": "2026-06-25T00:00:00Z",
      "ExpirationStr": "20260625",
      "BarcodeStr": "24062400390041,820-00002,2406240039,20260625"
    }
  },
  "LibraryType": "Linear",
  "RecipeValues": [
    {
      "Name": "filterMask",
      "Value": "R1:Y15N*-R2:Y15N*"
    }
  ],
  "AdvancedSettings": {
    "PolonyDensity": "HighDensity"
  }
}
""")
        with open(f"{run_path}/RunUploaded.json", "w") as stream:
            outcome = "OutcomeCompleted" if outcome_completed else "OutcomeFailed"
            stream.write(
                "{"
                + '"version":"1.0.0",'
                + '"instrument":"AV242106",'
                + '"instrumentId":"0000024023696901c5621014",'
                + '"runType":"Sequencing",'
                + '"recipeExecutionId":"rec.9590c80c95fc4eee8b3eb10c31251915",'
                + '"runID":"seq_66f5837f1ae1a35f10a2e594",'
                + f'"outcome":"{outcome}"'
                + "}"
            )

    if sync_finished:
        open(f"{run_path}/.sync_finished", "w").close()

    if demux_dir:
        os.mkdir(os.path.join(run_path, "Demultiplexing"))
        if demux_done:
            open(
                os.path.join(
                    run_path,
                    "Demultiplexing",
                    "RunStats.json",
                ),
                "w",
            ).close()
        for i in range(n_demux_subdirs):
            os.mkdir(os.path.join(run_path, f"Demultiplexing_{i}"))
            if demux_done:
                open(
                    os.path.join(
                        run_path,
                        f"Demultiplexing_{i}",
                        "RunStats.json",
                    ),
                    "w",
                ).close()

    return run_path


@mock.patch("taca.element.Element_Runs.ElementRunsConnection")
class TestRun:
    def test_init(self, mock_db: mock.Mock, create_dirs: pytest.fixture):
        tmp: tempfile.TemporaryDirectory = create_dirs
        run_dir = create_element_run_dir(tmp)

        run = to_test.Run(run_dir, get_config(tmp))
        assert run.run_dir == run_dir

    @pytest.mark.parametrize(
        "p",
        [
            {"run_finished": True, "outcome_completed": True, "expected": True},
            {"run_finished": True, "outcome_completed": False, "expected": False},
            {"run_finished": False, "outcome_completed": False, "expected": False},
        ],
        ids=["success", "failure", "ongoing"],
    )
    def test_check_sequencing_status(
        self,
        mock_db: mock.Mock,
        p: pytest.fixture,
        create_dirs: pytest.fixture,
    ):
        tmp: tempfile.TemporaryDirectory = create_dirs

        run = to_test.Run(
            create_element_run_dir(
                tmp,
                run_finished=p["run_finished"],
                outcome_completed=p["outcome_completed"],
            ),
            get_config(tmp),
        )
        assert run.check_sequencing_status() is p["expected"]

    @pytest.mark.parametrize(
        "p",
        [
            {"demux_dir": False, "demux_done": False, "expected": "not started"},
            {"demux_dir": True, "demux_done": False, "expected": "ongoing"},
            {"demux_dir": True, "demux_done": True, "expected": "finished"},
        ],
        ids=["not started", "ongoing", "finished"],
    )
    def test_get_demultiplexing_status(
        self, mock_db: mock.Mock, p: pytest.fixture, create_dirs: pytest.fixture
    ):
        tmp: tempfile.TemporaryDirectory = create_dirs

        run = to_test.Run(
            create_element_run_dir(
                tmp,
                demux_dir=p["demux_dir"],
                demux_done=p["demux_done"],
            ),
            get_config(tmp),
        )

        assert run.get_demultiplexing_status() == p["expected"]

    @pytest.mark.skip(reason="Not implemented yet")
    def test_generate_demux_command(self, mock_db):
        pass

    def test_start_demux(self, mock_db, create_dirs):
        tmp: tempfile.TemporaryDirectory = create_dirs
        with mock.patch("subprocess.Popen") as mock_Popen, mock.patch(
            "taca.element.Element_Runs.Run.generate_demux_command"
        ) as mock_command:
            mock_command.return_value = "test command"
            run = to_test.Run(create_element_run_dir(create_dirs), get_config(tmp))
            run.start_demux("mock_run_manifest", "mock_demux_dir")
            mock_command.assert_called_once_with("mock_run_manifest", "mock_demux_dir")
            mock_Popen.assert_called_once()

    @pytest.mark.skip(reason="Not implemented yet")
    def test_is_transferred(self, mock_db, create_dirs):
        pass

    @pytest.mark.skip(reason="Not implemented yet")
    def test_parse_rundir(self, mock_db, create_dirs):
        pass
