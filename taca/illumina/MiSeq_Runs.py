from taca.illumina.Standard_Runs import Standard_Runs


class MiSeq_Run(Standard_Runs):
    def __init__(self, run_dir, software, configuration):
        super(MiSeq_Run, self).__init__( run_dir, software, configuration)
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
