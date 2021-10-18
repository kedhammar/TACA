import os
import subprocess
import shutil
import glob
import logging

from taca.nanopore.nanopore import Nanopore
from taca.utils.config import CONFIG
from taca.utils.minion_barcodes import BARCODES

logger = logging.getLogger(__name__)

class MinION(Nanopore):
    """Minion run"""
    def __init__(self, run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet):
        super(MinION, self).__init__(run_dir)
        self.nanoseq_sample_sheet = nanoseq_sample_sheet
        self.anglerfish_sample_sheet = anglerfish_sample_sheet
        self.nanoseq_dir = os.path.join(self.run_dir, 'nanoseq_output')
        self.anglerfish_dir = os.path.join(self.run_dir, 'anglerfish_output')
        self.nanoseq_exit_status_file = os.path.join(self.run_dir, '.exitcode_for_nanoseq')
        self.anglerfish_exit_status_file = os.path.join(self.run_dir, '.exitcode_for_anglerfish')
        self.year_processed = self.run_id[0:4]
        self.flowcell_id = self.run_id.split('_')[3]

    def parse_lims_sample_sheet(self):
        """Generate nanoseq samplesheet based on Lims original."""
        self._get_original_samplesheet()
        if self.lims_samplesheet:
            self._parse_samplesheet()
        else:
            self.nanoseq_sample_sheet = ''
        return 

    def _get_original_samplesheet(self):
        """Find original lims sample sheet."""
        lims_samplesheet_dir = os.path.join(CONFIG.get('nanopore_analysis').get('samplesheets_dir'),
                                            self.year_processed)
        found_samplesheets = glob.glob(lims_samplesheet_dir + '/*' + self.flowcell_id + '*')
        if not found_samplesheets:
            logger.warn('No Lims sample sheets found for run {}'.format(self.run_id))
            self.lims_samplesheet = None
            return
        elif len(found_samplesheets) > 1:
            logger.warn('Found more than one Lims sample sheets for run {}'.format(self.run_id))
            self.lims_samplesheet = None
            return
        self.lims_samplesheet = found_samplesheets[0]
        return

    def _parse_samplesheet(self):
        """Parse Lims samplesheet into one suitable for nanoseq and anglerfish."""
        nanopore_kit = os.path.basename(self.lims_samplesheet).split('_')[0]
        self.nanoseq_sample_sheet = os.path.join(self.run_dir, nanopore_kit + '_sample_sheet.csv')
        self.anglerfish_sample_sheet = os.path.join(self.run_dir, 'anglerfish_sample_sheet.csv')
        nanoseq_content = 'sample,fastq,barcode,genome,transcriptome'
        anglerfish_content = ''
        with open(self.lims_samplesheet, 'r') as f:
            lines = sorted(f.readlines())
            first_sample_name = lines[0].split(',')[0]
            fastq_location = os.path.join(self.run_dir, 'nanoseq_output', 'guppy', 'fastq')
            for line in lines:
                sample_name, nanoseq_barcode, run_type, illumina_barcode = line.split(',')
                if nanoseq_barcode and nanoseq_barcode in BARCODES:
                    barcode = BARCODES[nanoseq_barcode]
                    is_pool = False
                else:
                    barcode = '0'
                    is_pool = True
                nanoseq_content += '\n' + sample_name + ',,' + barcode + ',,'  # Only need sample and barcode for now.
                if illumina_barcode:
                    # If there are no nanopore barcodes, the samples are from the same pool and will end up in
                    # the same nanoseq output file named after the firts sample in the sample sheet
                    if is_pool:
                        anglerfish_content += sample_name + ',' + run_type + ',' + illumina_barcode + ',' + os.path.join(fastq_location, first_sample_name + '.fastq.gz') + '\n'
                    # If the samples are not the same pool, the nanoseq output is named by the barcode
                    else:
                        anglerfish_content += sample_name + ',' + run_type + ',' + illumina_barcode + ',' + os.path.join(fastq_location, 'barcode' + barcode + '.fastq.gz') + '\n'
        with open(self.nanoseq_sample_sheet, 'w') as f:
            f.write(nanoseq_content)
        if anglerfish_content:
            with open(self.anglerfish_sample_sheet, 'w') as f:
                f.write(anglerfish_content)
        return

    def start_nanoseq(self):
        """Start Nanoseq analysis."""
        flowcell_product_code = self._get_flowcell_product_code() 
        kit_id = os.path.basename(self.nanoseq_sample_sheet).split('_')[0]
        if self._is_multiplexed():
            logger.info('Run {} is multiplexed. Starting nanoseq with --barcode_kit option'.format(self.run_dir))
            barcode_kit = self._get_barcode_kit()
            analysis_command = ('nextflow run nf-core/nanoseq --input ' + self.nanoseq_sample_sheet
                                + ' --input_path ' + os.path.join(self.run_dir, 'fast5')
                                + ' --outdir ' + os.path.join(self.run_dir, 'nanoseq_output')
                                + ' --flowcell ' + flowcell_product_code
                                + ' --guppy_gpu'
                                + ' --skip_alignment'
                                + ' --kit ' + kit_id
                                + ' --max_cpus 6'
                                + ' --max_memory 20.GB'
                                + ' --barcode_kit ' + barcode_kit
                                + ' -profile singularity; echo $? > .exitcode_for_nanoseq')
        else:
            logger.info('Run {} is not multiplexed. Starting nanoseq without --barcode_kit option'.format(self.run_dir))
            analysis_command = ('nextflow run nf-core/nanoseq --input ' + self.nanoseq_sample_sheet
                                + ' --input_path ' + os.path.join(self.run_dir, 'fast5')
                                + ' --outdir ' + os.path.join(self.run_dir, 'nanoseq_output')
                                + ' --flowcell ' + flowcell_product_code
                                + ' --guppy_gpu'
                                + ' --skip_alignment'
                                + ' --kit ' + kit_id
                                + ' --max_cpus 6'
                                + ' --max_memory 20.GB'
                                + ' -profile singularity; echo $? > .exitcode_for_nanoseq')

        try:
            p_handle = subprocess.Popen(analysis_command, stdout=subprocess.PIPE, shell=True, cwd=self.run_dir)
            logger.info('Started Nanoseq for run {}'.format(self.run_dir))
        except subprocess.CalledProcessError:
            logger.warn('An error occurred while starting the Nanoseq for run {}. '
                        'Please check the logfile for info.'.format(self.run_dir))
        return

    def _get_flowcell_product_code(self):
        """Look for flow_cell_product_code in report.md and return the corresponding value."""
        report_file = glob.glob(self.run_dir + '/report*.md')[0]
        with open(report_file, 'r') as f:
            for line in f.readlines():
                if 'flow_cell_product_code' in line:
                    return line.split('"')[3]

    def _is_multiplexed(self):
        """Look in the sample_sheet and return True if the run was multiplexed, else False.
        Assumes that a run that has not been multiplexed has the barcode 0."""
        with open(self.nanoseq_sample_sheet, 'r') as f:
            for i, line in enumerate(f):
                if i == 1:  # Only need to check first non-header line
                    line_entries = line.split(',')
        if line_entries[2] == '0':
            return False
        else:
            return True

    def _get_barcode_kit(self):
        """Figure out which barcode kit was used. Assumes only one kit is ever used."""
        with open(self.nanoseq_sample_sheet, 'r') as f:
            for i, line in enumerate(f):
                if i == 1:  # Only need to check first non-header line
                    line_entries = line.split(',')
        if int(line_entries[2]) <= 12:
            return 'EXP-NBD104'
        elif int(line_entries[2]) >= 13:
            return 'EXP-NBD114'

    def check_exit_status(self, status_file):
        """Read pipeline exit status file and return True if 0, False if anything else"""
        with open(status_file, 'r') as f:
            exit_status = f.readline().strip()
        return exit_status == '0'

    def start_anglerfish(self):
        """Start Anglerfish."""
        os.makedirs(self.anglerfish_dir)
        anglerfish_command = ('anglerfish.py'
                            + ' --samplesheet ' + self.anglerfish_sample_sheet
                            + ' --out_fastq ' + self.anglerfish_dir
                            + ' --threads 2'
                            + ' --skip_demux'
                            + ' --skip_fastqc; echo $? > .exitcode_for_anglerfish')
        try:
            p_handle = subprocess.Popen(anglerfish_command, stdout=subprocess.PIPE, shell=True, cwd=self.run_dir)
            logger.info('Started Anglerfish for run {} using: {}'.format(self.run_dir, anglerfish_command))
        except subprocess.CalledProcessError:
            logger.warn('An error occurred while starting the Anglerfish for run {}. '
                        'Please check the logfile for info.'.format(self.run_dir))
        return

    def copy_results_for_lims(self):
        """Find results and copy to lims directory."""
        year_processed = self.run_id[0:4]
        lims_result_file = os.path.join(CONFIG.get('nanopore_analysis').get('lims_results_dir'),
                                        year_processed, 'anglerfish_stats_' + self.flowcell_id + '.txt')
        anglerfish_results = self._find_anglerfish_results()
        try:
            shutil.copyfile(anglerfish_results, lims_result_file)
        except OSError as e:
            logger.warn('An error occurred while copying the Anglerfish results for {} to lims: {}'.format(self.run_id, e))
        return

    def _find_anglerfish_results(self):
        """Return location of Anglerfish results."""
        results_file = ''
        for sub_dir in os.listdir(self.anglerfish_dir):
            if 'anglerfish_stats.txt' in os.listdir(os.path.join(self.anglerfish_dir, sub_dir)):
                results_file = os.path.join(self.anglerfish_dir, sub_dir, 'anglerfish_stats.txt')
                return results_file
        if not results_file:
            logger.warn('Could not find any Anglerfish results in {}'.format(self.anglerfish_dir))
