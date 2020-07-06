import os
from taca.illumina.HiSeq_Runs import HiSeq_Run
from flowcell_parser.classes import SampleSheetParser


class MiSeq_Run(HiSeq_Run):

    def __init__(self,  path_to_run, configuration):
        #constructor, it returns a MiSeq object only if the MiSeq run belongs to NGI facility, i.e., contains
        #Application or production in the Description
        super(MiSeq_Run, self).__init__( path_to_run, configuration)
        self._set_sequencer_type()
        self._set_run_type()
        self._copy_samplesheet()

    def _set_sequencer_type(self):
        self.sequencer_type = "MiSeq"

    def _set_run_type(self):
        ssname = os.path.join(self.run_dir, 'Data', 'Intensities', 'BaseCalls','SampleSheet.csv')
        if not os.path.exists(ssname):
            #case in which no samplesheet is found, assume it is a non NGI run
            self.run_type = "NON-NGI-RUN"
        else:
            #it SampleSheet exists try to see if it is a NGI-run
            ssparser = SampleSheetParser(ssname)
            if ssparser.header['Description'] == "Production" or ssparser.header['Description'] == "Applications":
                self.run_type = "NGI-RUN"
            else:
            #otherwise this is a non NGI run
                self.run_type = "NON-NGI-RUN"

    def _get_samplesheet(self):
        """
        Locate and parse the samplesheet for a run.
        In MiSeq case this is located in FC_DIR/Data/Intensities/BaseCalls/SampleSheet.csv
        """
        ssname = os.path.join(self.run_dir, 'Data', 'Intensities', 'BaseCalls','SampleSheet.csv')
        if os.path.exists(ssname):
            #if exists parse the SampleSheet
            return ssname
        else:
            #some MiSeq runs do not have the SampleSheet at all, in this case assume they are non NGI.
            #not real clean solution but what else can be done if no samplesheet is provided?
            return None

    def _generate_clean_samplesheet(self, ssparser):
        """
        Will generate a 'clean' samplesheet, for bcl2fastq
        """
        output=""
        #Header
        output+="[Header]{}".format(os.linesep)
        for field in ssparser.header:
            output+="{},{}".format(field.rstrip(), ssparser.header[field].rstrip())
            output+=os.linesep
        #now parse the data section
        data = []
        for line in ssparser.data:
            entry = {}
            for field, value in line.items():
                if ssparser.dfield_sid in field:
                    entry[field] ='Sample_{}'.format(value)
                elif ssparser.dfield_proj in field:
                    entry[field] = value.replace(".", "__")
                else:
                    entry[field] = value
            if 'Lane' not in entry:
                entry['Lane'] = '1'
            if 'index' not in entry:
                entry['index'] = ''
            if 'I7_Index_ID' not in entry:
                entry['I7_Index_ID'] = ''
            if 'index2' not in entry:
                entry['index2'] = ''
            if 'I5_Index_ID' not in entry:
                entry['I5_Index_ID'] = ''
            data.append(entry)

        fields_to_output = ['Lane', ssparser.dfield_sid, ssparser.dfield_snm, 'index', ssparser.dfield_proj, 'I7_Index_ID', 'index2', 'I5_Index_ID']
        #now create the new SampleSheet data section
        output+="[Data]{}".format(os.linesep)
        for field in ssparser.datafields:
            if field not in fields_to_output:
                fields_to_output.append(field)
        output+=",".join(fields_to_output)
        output+=os.linesep
        #now process each data entry and output it
        for entry in data:
            line = []
            for field in fields_to_output:
                line.append(entry[field])
            output+=",".join(line)
            output+=os.linesep
        return output
