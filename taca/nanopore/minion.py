from taca.nanopore.nanopore import Nanopore

class MinION(Nanopore):
    """Minion run"""
    
    def __init__(self):
        pass

    def parse_lims_sample_sheet(run_dir):
        """Generate nanoseq samplesheet based on Lims original."""
        run_id = os.path.basename(run_dir)
        lims_samplesheet = get_original_samplesheet(run_id)
        if lims_samplesheet:
            nanoseq_samplesheet_location = parse_samplesheet(run_dir, lims_samplesheet)
        else:
            nanoseq_samplesheet_location = ''
            anglerfish_samplesheet_location = ''
        return nanoseq_samplesheet_location

    def get_original_samplesheet(run_id):
        """Find original lims sample sheet."""
        year_processed = run_id[0:4]
        flowcell_id = run_id.split('_')[3]
        lims_samplesheet_dir = os.path.join(CONFIG.get('nanopore_analysis').get('samplesheets_dir'),
                                            year_processed)
        found_samplesheets = glob.glob(lims_samplesheet_dir + '/*'+ flowcell_id + '*')
        if not found_samplesheets:
            logger.warn('No Lims sample sheets found for run {}'.format(run_id))
            return
        elif len(found_samplesheets) > 1:
            logger.warn('Found more than one Lims sample sheets for run {}'.format(run_id))
            return
        return found_samplesheets[0]

    def parse_samplesheet(run_dir, lims_samplesheet):
        """Parse Lims samplesheet into one suitable for nanoseq and anglerfish."""
        nanopore_kit = os.path.basename(lims_samplesheet).split('_')[0]
        nanoseq_samplesheet = os.path.join(run_dir, nanopore_kit + '_sample_sheet.csv')
        anglerfish_samplesheet = os.path.join(run_dir, 'anglerfish_sample_sheet.csv')
        nanoseq_content = 'sample,fastq,barcode,genome,transcriptome'
        anglerfish_content = ''
        with open(lims_samplesheet, 'r') as f:
            lines = sorted(f.readlines())
            first_sample_name = lines[0].split(',')[0] # Get name of first sample
            fastq_location = os.path.join(run_dir, 'nanoseq_output', 'guppy', 'fastq') # Set the location of the first sample
            for line in lines:
                sample_name, nanoseq_barcode, run_type, illumina_barcode, location = line.split(',') #TODO remove location once/if removed in lims
                if nanoseq_barcode and nanoseq_barcode in BARCODES:
                    barcode = BARCODES[nanoseq_barcode]
                    is_pool = False
                else:
                    barcode = '0'
                    is_pool = True
                nanoseq_content += '\n' + sample_name + ',,' + barcode + ',,' # Only need sample and barcode for now.
                if illumina_barcode:
                    # If there are no nanopore barcodes, the samples are from the same pool and will end up in
                    # the same nanoseq output file named after the firts sample in the sample sheet
                    if is_pool:
                        anglerfish_content += sample_name + ',' + run_type + ',' + illumina_barcode + ',' + os.path.join(fastq_location, first_sample_name + '.fastq.gz') + '\n'
                    # If the samples are not the same pool, the nanoseq output is named by the barcode
                    else:
                        anglerfish_content += sample_name + ',' + run_type + ',' + illumina_barcode + ',' + os.path.join(fastq_location, 'barcode' + barcode + '.fastq.gz') + '\n'
        with open(nanoseq_samplesheet, 'w') as f:
            f.write(nanoseq_content)
        if anglerfish_content:
            with open(anglerfish_samplesheet, 'w') as f:
                f.write(anglerfish_content)
        return nanoseq_samplesheet

    def start_nanoseq(run_dir, sample_sheet):
        """Start Nanoseq analysis."""
        flowcell_id = get_flowcell_id(run_dir)
        kit_id = os.path.basename(sample_sheet).split('_')[0]
        if is_multiplexed(sample_sheet):
            logger.info('Run {} is multiplexed. Starting nanoseq with --barcode_kit option'.format(run_dir))
            barcode_kit = get_barcode_kit(sample_sheet)
            analysis_command = ('nextflow run nf-core/nanoseq --input ' + sample_sheet
                                + ' --input_path ' + os.path.join(run_dir, 'fast5')
                                + ' --outdir ' + os.path.join(run_dir, 'nanoseq_output')
                                + ' --flowcell ' + flowcell_id
                                + ' --guppy_gpu'
                                + ' --skip_alignment'
                                + ' --kit ' + kit_id
                                + ' --max_cpus 6'
                                + ' --max_memory 20.GB'
                                + ' --barcode_kit ' + barcode_kit
                                + ' -profile singularity; echo $? > .exitcode_for_nanoseq')
        else:
            logger.info('Run {} is not multiplexed. Starting nanoseq without --barcode_kit option'.format(run_dir))
            analysis_command = ('nextflow run nf-core/nanoseq --input ' + sample_sheet
                                + ' --input_path ' + os.path.join(run_dir, 'fast5')
                                + ' --outdir ' + os.path.join(run_dir, 'nanoseq_output')
                                + ' --flowcell ' + flowcell_id
                                + ' --guppy_gpu'
                                + ' --skip_alignment'
                                + ' --kit ' + kit_id
                                + ' --max_cpus 6'
                                + ' --max_memory 20.GB'
                                + ' -profile singularity; echo $? > .exitcode_for_nanoseq')

        try:
            p_handle = subprocess.Popen(analysis_command, stdout=subprocess.PIPE, shell=True, cwd=run_dir)
            logger.info('Started Nanoseq for run {}'.format(run_dir))
        except subprocess.CalledProcessError:
            logger.warn('An error occurred while starting the Nanoseq for run {}. '
                        'Please check the logfile for info.'.format(run_dir))
        return

    def get_flowcell_id(run_dir):
        """Look for flow_cell_product_code in report.md and return the corresponding value."""
        report_file = glob.glob(run_dir + '/report*.md')[0]
        with open(report_file, 'r') as f:
            for line in f.readlines():
                if 'flow_cell_product_code' in line:
                    return line.split('"')[3]

    def is_multiplexed(sample_sheet):
        """Look in the sample_sheet and return True if the run was multiplexed, else False.
        Assumes that a run that has not been multiplexed has the barcode 0."""
        with open(sample_sheet, 'r') as f:
            for i, line in enumerate(f):
                if i == 1: # Only need to check first non-header line
                    line_entries = line.split(',')
        if line_entries[2] == '0':
            return False
        else:
            return True

    def get_barcode_kit(sample_sheet):
        """Figure out which barcode kit was used. Assumes only one kit is ever used."""
        with open(sample_sheet, 'r') as f:
            for i, line in enumerate(f):
                if i == 1: # Only need to check first non-header line
                    line_entries = line.split(',')
        if int(line_entries[2]) <= 12:
            return 'EXP-NBD104'
        elif int(line_entries[2]) >= 13:
            return 'EXP-NBD114'
        barcode_kit = get_barcode_kit(sample_sheet)

    def check_exit_status(status_file):
        """Read pipeline exit status file and return True if 0, False if anything else"""
        with open(status_file, 'r') as f:
            exit_status = f.readline().strip()
        return exit_status == '0'

    def start_anglerfish(run_dir, af_sample_sheet, output_dir):
        """Start Anglerfish."""
        os.makedirs(output_dir)
        anglerfish_command = ('anglerfish.py'
                            + ' --samplesheet ' + af_sample_sheet
                            + ' --out_fastq ' + output_dir
                            + ' --threads 2'
                            + ' --skip_demux'
                            + ' --skip_fastqc; echo $? > .exitcode_for_anglerfish')
        try:
            p_handle = subprocess.Popen(anglerfish_command, stdout=subprocess.PIPE, shell=True, cwd=run_dir)
            logger.info('Started Anglerfish for run {} using: {}'.format(run_dir, anglerfish_command))
        except subprocess.CalledProcessError:
            logger.warn('An error occurred while starting the Anglerfish for run {}. '
                        'Please check the logfile for info.'.format(run_dir))
        return

    def copy_results_for_lims(run_dir, anglerfish_results_dir):
        """Find results and copy to lims directory."""
        run_id = os.path.basename(run_dir)
        year_processed = run_id[0:4]
        flowcell_id = run_id.split('_')[3]
        lims_result_file = os.path.join(CONFIG.get('nanopore_analysis').get('lims_results_dir'),
                                        year_processed, 'anglerfish_stats_' + flowcell_id + '.txt')
        anglerfish_results = find_anglerfish_results(anglerfish_results_dir)
        try:
            shutil.copyfile(anglerfish_results, lims_result_file)
        except OSError as e:
            logger.warn('An error occurred while copying the Anglerfish results for {} to lims: {}'.format(run_id, e))
        return

    def find_anglerfish_results(anglerfish_dir):
        """Return location of Anglerfish results."""
        results_file = ''
        for sub_dir in os.listdir(anglerfish_dir):
            if 'anglerfish_stats.txt' in os.listdir(os.path.join(anglerfish_dir, sub_dir)):
                results_file = os.path.join(anglerfish_dir, sub_dir, 'anglerfish_stats.txt')
                return results_file
        if not results_file:
            logger.warn('Could not find any Anglerfish results in {}'.format(anglerfish_dir))
