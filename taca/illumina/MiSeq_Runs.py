import os
import shutil
import logging
from taca.illumina.HiSeq_Runs import HiSeq_Run
from flowcell_parser.classes import SampleSheetParser

logger = logging.getLogger(__name__)

class MiSeq_Run(HiSeq_Run):

    def __init__(self,  path_to_run, configuration):
        # Constructor, it returns a MiSeq object only if the MiSeq run belongs to NGI facility, i.e., contains
        # Application or production in the Description
        super(MiSeq_Run, self).__init__( path_to_run, configuration)
        self._set_sequencer_type()
        self._set_run_type()
        self._copy_samplesheet()

    def _set_sequencer_type(self):
        self.sequencer_type = 'MiSeq'

    def _set_run_type(self):
        ssname = os.path.join(self.run_dir,
                              'SampleSheet.csv')
        if not os.path.exists(ssname):
            # Case in which no samplesheet is found, assume it is a non NGI run
            self.run_type = 'NON-NGI-RUN'
        else:
            # If SampleSheet exists try to see if it is a NGI-run
            ssparser = SampleSheetParser(ssname)
            if ssparser.header['Description'] == 'Production' or ssparser.header['Description'] == 'Applications':
                self.run_type = 'NGI-RUN'
            else:
                # Otherwise this is a non NGI run
                self.run_type = 'NON-NGI-RUN'

    def _get_samplesheet(self):
        """Locate and parse the samplesheet for a run.
        In MiSeq case this is located in FC_DIR/SampleSheet.csv
        """
        ssname = os.path.join(self.run_dir,
                              'SampleSheet.csv')
        if os.path.exists(ssname):
            # If exists parse the SampleSheet
            return ssname
        else:
            # Some MiSeq runs do not have the SampleSheet at all, in this case assume they are non NGI.
            # Not real clean solution but what else can be done if no samplesheet is provided?
            return None

    def _copy_samplesheet(self):
        ssname = self._get_samplesheet()
        if ssname is None:
            return None
        ssparser = SampleSheetParser(ssname)
        # Copy the original samplesheet locally.
        # Copy again if already done as there might have been changes to the samplesheet
        try:
            shutil.copy(ssname, os.path.join(self.run_dir, '{}.csv'.format(self.flowcell_id)))
            ssname = os.path.join(self.run_dir, os.path.split(ssname)[1])
        except:
            raise RuntimeError('unable to copy file {} to destination {}'.format(ssname, self.run_dir))

        # This sample sheet has been created by the LIMS and copied by a sequencing operator. It is not ready
        # to be used it needs some editing.
        # This will contain the samplesheet with all the renaiming to be used with bcl2fastq
        samplesheet_dest = os.path.join(self.run_dir, 'SampleSheet_copy.csv')
        # Check that the samplesheet is not already present. In this case go the next step
        if os.path.exists(samplesheet_dest):
            logger.info('SampleSheet_copy.csv found ... overwriting it')
        try:
            with open(samplesheet_dest, 'w') as fcd:
                fcd.write(self._generate_clean_samplesheet(ssparser))
        except Exception as e:
            logger.error(e)
            return False
        logger.info(('Created SampleSheet_copy.csv for Flowcell {} in {} '.format(self.id, samplesheet_dest)))
        # SampleSheet.csv generated
        # When demultiplexing SampleSheet.csv is the one I need to use
        self.runParserObj.samplesheet  = SampleSheetParser(os.path.join(self.run_dir, 'SampleSheet_copy.csv'))
        if not self.runParserObj.obj.get('samplesheet_csv'):
            self.runParserObj.obj['samplesheet_csv'] = self.runParserObj.samplesheet.data

    def _generate_clean_samplesheet(self, ssparser):
        """Will generate a 'clean' samplesheet, for bcl2fastq"""
        output = u''
        # Header
        output += '[Header]{}'.format(os.linesep)
        for field in sorted(ssparser.header):
            output += '{},{}'.format(field.rstrip(), ssparser.header[field].rstrip())
            output += os.linesep
        # Parse the data section
        data = []
        for line in ssparser.data:
            entry = {}
            for field, value in line.items():
                if ssparser.dfield_sid in field:
                    entry[field] = 'Sample_{}'.format(value)
                elif ssparser.dfield_proj in field:
                    entry[field] = value.replace('.', '__')
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

        fields_to_output = ['Lane',
                            ssparser.dfield_sid,
                            ssparser.dfield_snm,
                            'index',
                            ssparser.dfield_proj,
                            'I7_Index_ID',
                            'index2',
                            'I5_Index_ID']
        # Create the new SampleSheet data section
        output += '[Data]{}'.format(os.linesep)
        for field in ssparser.datafields:
            if field not in fields_to_output:
                fields_to_output.append(field)
        output += ','.join(fields_to_output)
        output += os.linesep
        # Process each data entry and output it
        for entry in data:
            line = []
            for field in fields_to_output:
                line.append(entry[field])
            output += ','.join(line)
            output += os.linesep
        return output
