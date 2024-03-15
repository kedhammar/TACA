import logging
import os
import re
import shutil

from flowcell_parser.classes import SampleSheetParser

from taca.illumina.Standard_Runs import Standard_Run

logger = logging.getLogger(__name__)

TENX_SINGLE_PAT = re.compile("SI-(?:GA|NA)-[A-H][1-9][0-2]?")
TENX_DUAL_PAT = re.compile("SI-(?:TT|NT|NN|TN|TS)-[A-H][1-9][0-2]?")
SMARTSEQ_PAT = re.compile("SMARTSEQ[1-9]?-[1-9][0-9]?[A-P]")
IDT_UMI_PAT = re.compile("([ATCG]{4,}N+$)")
RECIPE_PAT = re.compile("[0-9]+-[0-9]+")


class MiSeq_Run(Standard_Run):
    def __init__(self, run_dir, software, configuration):
        super().__init__(run_dir, software, configuration)
        self._set_sequencer_type()
        self._set_run_type()
        self._copy_samplesheet()

    def _set_sequencer_type(self):
        self.sequencer_type = "MiSeq"

    def _set_run_type(self):
        self.run_type = "NGI-RUN"

    def _get_samplesheet(self):
        """Locate and parse the samplesheet for a run.
        In MiSeq case this is located in FC_DIR/SampleSheet.csv
        """
        ssname = os.path.join(self.run_dir, "SampleSheet.csv")
        if os.path.exists(ssname):
            # If exists parse the SampleSheet
            return ssname
        else:
            # Some MiSeq runs do not have the SampleSheet at all, in this case assume they are non NGI.
            # Not real clean solution but what else can be done if no samplesheet is provided?
            return None

    def _copy_samplesheet(self):
        ssname = self._get_samplesheet()
        runSetup = self.runParserObj.runinfo.get_read_configuration()
        # Load index files
        indexfile = dict()
        try:
            indexfile["tenX"] = self.CONFIG[self.software]["tenX_index_path"]
        except KeyError:
            logger.error("Path to index file (10X) not found in the config file")
            raise RuntimeError
        try:
            indexfile["smartseq"] = self.CONFIG[self.software]["smartseq_index_path"]
        except KeyError:
            logger.error("Path to index file (Smart-seq) not found in the config file")
            raise RuntimeError
        if ssname is None:
            return None
        ssparser = SampleSheetParser(ssname)
        self.sample_table = self._classify_samples(indexfile, ssparser, runSetup)
        # Copy the original samplesheet locally.
        # Copy again if already done as there might have been changes to the samplesheet
        try:
            shutil.copy(ssname, os.path.join(self.run_dir, f"{self.flowcell_id}.csv"))
            ssname = os.path.join(self.run_dir, os.path.split(ssname)[1])
        except:
            raise RuntimeError(
                f"unable to copy file {ssname} to destination {self.run_dir}"
            )

        # This sample sheet has been created by the LIMS and copied by a sequencing operator. It is not ready
        # to be used it needs some editing.
        # This will contain the samplesheet with all the renaiming to be used with bcl2fastq
        samplesheet_dest = os.path.join(self.run_dir, "SampleSheet_copy.csv")
        # Check that the samplesheet is not already present. In this case go the next step
        if os.path.exists(samplesheet_dest):
            logger.info("SampleSheet_copy.csv found ... overwriting it")
        try:
            with open(samplesheet_dest, "w") as fcd:
                fcd.write(
                    self._generate_clean_samplesheet(
                        ssparser,
                        indexfile,
                        fields_to_remove=None,
                        rename_samples=True,
                        rename_qPCR_suffix=True,
                        fields_qPCR=[ssparser.dfield_snm],
                    )
                )
        except Exception as e:
            logger.error(e)
            return False
        logger.info(
            f"Created SampleSheet_copy.csv for Flowcell {self.id} in {samplesheet_dest} "
        )
        # SampleSheet.csv generated
        # When demultiplexing SampleSheet.csv is the one I need to use
        self.runParserObj.samplesheet = SampleSheetParser(
            os.path.join(self.run_dir, "SampleSheet_copy.csv")
        )
        if not self.runParserObj.obj.get("samplesheet_csv"):
            self.runParserObj.obj["samplesheet_csv"] = (
                self.runParserObj.samplesheet.data
            )

    def _generate_clean_samplesheet(
        self,
        ssparser,
        indexfile,
        fields_to_remove=None,
        rename_samples=True,
        rename_qPCR_suffix=False,
        fields_qPCR=None,
    ):
        """Generate a 'clean' samplesheet, the given fields will be removed.
        If rename_samples is True, samples prepended with 'Sample_'  are renamed to match the sample name
        Will also replace 10X or Smart-seq indicies (e.g. SI-GA-A3 into TGTGCGGG)
        Note that the index 2 of 10X or Smart-seq dual indexes will be converted to RC
        """
        output = ""
        compl = {"A": "T", "C": "G", "G": "C", "T": "A"}
        # Expand the ssparser if there are lanes with 10X or Smart-seq samples
        index_dict_tenX = self._parse_10X_indexes(indexfile["tenX"])
        index_dict_smartseq = self._parse_smartseq_indexes(indexfile["smartseq"])
        # Replace 10X or Smart-seq indices
        for sample in ssparser.data:
            if sample["index"] in index_dict_tenX.keys():
                tenX_index = sample["index"]
                # In the case of 10X dual indexes, replace index and index2
                if TENX_DUAL_PAT.findall(tenX_index):
                    sample["index"] = index_dict_tenX[tenX_index][0]
                    sample["index2"] = "".join(
                        reversed(
                            [
                                compl.get(b, b)
                                for b in index_dict_tenX[tenX_index][1]
                                .replace(",", "")
                                .upper()
                            ]
                        )
                    )
                # In the case of 10X single indexes, replace the index name with the 4 actual indicies
                else:
                    x = 0
                    indices_number = len(index_dict_tenX[tenX_index])
                    while x < indices_number - 1:
                        new_sample = dict(sample)
                        new_sample["index"] = index_dict_tenX[tenX_index][x]
                        ssparser.data.append(new_sample)
                        x += 1
                    # Set the original 10X index to the 4th correct index
                    sample["index"] = index_dict_tenX[tenX_index][x]
            elif SMARTSEQ_PAT.findall(sample["index"]):
                x = 0
                smartseq_index = sample["index"].split("-")[1]
                indices_number = len(index_dict_smartseq[smartseq_index])
                while x < indices_number - 1:
                    new_sample = dict(sample)
                    new_sample["index"] = index_dict_smartseq[smartseq_index][x][0]
                    new_sample["index2"] = "".join(
                        reversed(
                            [
                                compl.get(b, b)
                                for b in index_dict_smartseq[smartseq_index][x][1]
                                .replace(",", "")
                                .upper()
                            ]
                        )
                    )
                    ssparser.data.append(new_sample)
                    x += 1
                sample["index"] = index_dict_smartseq[smartseq_index][x][0]
                sample["index2"] = "".join(
                    reversed(
                        [
                            compl.get(b, b)
                            for b in index_dict_smartseq[smartseq_index][x][1]
                            .replace(",", "")
                            .upper()
                        ]
                    )
                )

        # Sort to get the added indicies from 10x in the right place
        # Python 3 doesn't support sorting a list of dicts implicitly. Sort by lane and then Sample_ID
        ssparser.data.sort(key=lambda item: (item.get("Lane"), item.get("Sample_ID")))

        if not fields_to_remove:
            fields_to_remove = []
        # Header
        output += f"[Header]{os.linesep}"
        for field in sorted(ssparser.header):
            output += f"{field.rstrip()},{ssparser.header[field].rstrip()}"
            output += os.linesep
        # Data
        output += f"[Data]{os.linesep}"
        datafields = []
        for field in ssparser.datafields:
            if field not in fields_to_remove:
                datafields.append(field)
        output += ",".join(datafields)
        output += os.linesep
        for line in ssparser.data:
            line_ar = []
            for field in datafields:
                value = line[field]
                if rename_samples and ssparser.dfield_sid in field:
                    try:
                        if rename_qPCR_suffix and ssparser.dfield_snm in fields_qPCR:
                            # Substitute SampleID with SampleName, add Sample_ as prefix and remove __qPCR_ suffix
                            value = re.sub(
                                "__qPCR_$", "", f"Sample_{line[ssparser.dfield_snm]}"
                            )
                        else:
                            # Substitute SampleID with SampleName, add Sample_ as prefix
                            value = f"Sample_{line[ssparser.dfield_snm]}"
                    except:
                        # Otherwise add Sample_ as prefix
                        value = f"Sample_{line[ssparser.dfield_sid]}"
                elif rename_qPCR_suffix and field in fields_qPCR:
                    value = re.sub("__qPCR_$", "", line[field])
                line_ar.append(value)
            output += ",".join(line_ar)
            output += os.linesep
        return output
