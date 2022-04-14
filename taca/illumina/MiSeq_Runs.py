import os
import re
import shutil
import logging
from taca.illumina.HiSeq_Runs import HiSeq_Run
from flowcell_parser.classes import SampleSheetParser

logger = logging.getLogger(__name__)

IDX_BM_PAT = re.compile('I[0-9]*')

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


    def _generate_bcl2fastq_command(self, base_masks, strict=True, suffix=0, mask_short_adapter_reads=False):
        """Generates the command to demultiplex with the given base_masks.
        If strict is set to true demultiplex only lanes in base_masks
        """
        logger.info('Building bcl2fastq command')
        cl = [self.CONFIG.get('bcl2fastq')['bin']]
        if 'options' in self.CONFIG.get('bcl2fastq'):
            cl_options = self.CONFIG['bcl2fastq']['options']
            # Append all options that appear in the configuration file to the main command.
            for option in cl_options:
                if isinstance(option, dict):
                    opt, val = list(option.items())[0]
                    # Skip output-dir has I might need more than one
                    if 'output-dir' not in opt:
                        cl.extend(['--{}'.format(opt), str(val)])
                else:
                    cl.append('--{}'.format(option))

        # Add the base_mask for each lane
        tiles = []
        samplesheetMaskSpecific = os.path.join(os.path.join(self.run_dir, 'SampleSheet_{}.csv'.format(suffix)))
        output_dir = 'Demultiplexing_{}'.format(suffix)
        cl.extend(['--output-dir', output_dir])

        with open(samplesheetMaskSpecific, 'w') as ssms:
            ssms.write(u'[Header]\n')
            ssms.write(u'[Data]\n')
            ssms.write(u','.join(self.runParserObj.samplesheet.datafields))
            ssms.write(u'\n')
            for lane in sorted(base_masks):
                # Iterate thorugh each lane and add the correct --use-bases-mask for that lane
                # There is a single basemask for each lane, I checked it a couple of lines above
                base_mask = [base_masks[lane][bm]['base_mask'] for bm in base_masks[lane]][0] # Get the base_mask
                # Add the extra command option if we have samples with single short index
                idx_bm = [x[0] for x in [IDX_BM_PAT.findall(bm) for bm in base_mask] if len(x)>0]
                if len(idx_bm)==1:
                    if int(re.findall('\d+',idx_bm[0])[0]) <= 6:
                         cl_options.extend(self.CONFIG['bcl2fastq']['options_short_single_index'])
                base_mask_expr = '{}:'.format(lane) + ','.join(base_mask)
                cl.extend(['--use-bases-mask', base_mask_expr])
                if strict:
                    tiles.extend(['s_{}'.format(lane)])
                # These are all the samples that need to be demux with this samplemask in this lane
                samples   = [base_masks[lane][bm]['data'] for bm in base_masks[lane]][0]
                for sample in samples:
                    for field in self.runParserObj.samplesheet.datafields:
                        if field == 'index' and 'NOINDEX' in sample[field].upper():
                            ssms.write(u',') # This is emtpy due to NoIndex issue
                        else:
                            ssms.write(u'{},'.format(sample[field]))
                    ssms.write(u'\n')
            if strict:
                cl.extend(['--tiles', ','.join(tiles)])
        cl.extend(['--sample-sheet', samplesheetMaskSpecific])
        if mask_short_adapter_reads:
            cl.extend(['--mask-short-adapter-reads', '0'])

        logger.info(('BCL to FASTQ command built {} '.format(' '.join(cl))))
        return cl
