import os
import subprocess
import shutil
import glob
import logging
import pathlib

from taca.nanopore.nanopore import Nanopore
from taca.utils.config import CONFIG
from taca.utils.minion_barcodes import BARCODES

logger = logging.getLogger(__name__)

class MinIONqc(Nanopore):
    """Minion QC run"""
    def __init__(self, run_dir, nanoseq_sample_sheet, anglerfish_sample_sheet):
        super(MinIONqc, self).__init__(run_dir)
        self.transfer_details = CONFIG.get('nanopore_analysis').get('minion_qc_run').get('transfer')
        self.transfer_log = self.transfer_details.get('transfer_file')
        self.archive_dir = self.transfer_details.get('finished_dir')
        self.anglerfish_sample_sheet = anglerfish_sample_sheet
               
        self.basecalled_dir = os.path.join(self.run_dir, 'restuls?') #TODO: get actual name of results dir
        self.anglerfish_dir = os.path.join(self.run_dir, 'anglerfish_output')
        self.nanoseq_exit_status_file = os.path.join(self.run_dir, '.exitcode_for_nanoseq') #TODO: replace with run finished indicator if needed
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

    def _get_original_samplesheet(self):
        """Find original lims sample sheet."""
        lims_samplesheet_dir = os.path.join(CONFIG.get('nanopore_analysis').get('minion_qc_run').get('samplesheets_dir'),
                                            self.year_processed)
        found_samplesheets = glob.glob(lims_samplesheet_dir + '/*' + self.flowcell_id + '*')
        if not found_samplesheets:
            logger.warn('No Lims sample sheets found for run {}. Skipping it.'.format(self.run_id))
            self.lims_samplesheet = None
        elif len(found_samplesheets) > 1:
            logger.warn('Found more than one Lims sample sheets for run {}. Skipping it.'.format(self.run_id))
            self.lims_samplesheet = None
        else:
            self.lims_samplesheet = found_samplesheets[0]

    def _parse_samplesheet(self):
        """Parse Lims samplesheet into one suitable for nanoseq and anglerfish."""
        nanopore_kit = os.path.basename(self.lims_samplesheet).split('_')[1]
        self.nanoseq_sample_sheet = os.path.join(self.run_dir, nanopore_kit + '_sample_sheet.csv')
        self.anglerfish_sample_sheet = os.path.join(self.run_dir, 'anglerfish_sample_sheet.csv')
        nanoseq_content = 'group,replicate,barcode,input_file,genome,transcriptome' #'sample,fastq,barcode,genome,transcriptome'
        anglerfish_content = ''
        with open(self.lims_samplesheet, 'r') as f:
            lines = sorted(f.readlines())
            first_sample_name = lines[0].split(',')[0]
            fastq_location = os.path.join(self.run_dir, 'nanoseq_output', 'guppy', 'fastq')
            pool_barcodes = []
            for line in lines:
                sample_name, nanoseq_barcode, run_type, illumina_barcode = line.split(',')
                illumina_barcode = illumina_barcode.strip('\n')
                if nanoseq_barcode and nanoseq_barcode in BARCODES:
                    barcode = BARCODES[nanoseq_barcode]
                    is_single_pool = False
                else:
                    barcode = '0'
                    is_single_pool = True

                if barcode not in pool_barcodes:  # If there are multiple pools they should be treated like one sample each in nanoseq
                    nanoseq_content += '\n' + sample_name + ',1,' + barcode + ',,,'
                    pool_barcodes.append(barcode)

                if illumina_barcode:
                    # If there are no nanopore barcodes, the samples are from the same pool and will end up in
                    # the same nanoseq output file named after the first sample in the sample sheet
                    if is_single_pool:
                        anglerfish_content += sample_name + ',' + run_type + ',' + illumina_barcode + ',' + os.path.join(fastq_location, first_sample_name + '.fastq.gz') + '\n'
                    # If the samples are not the same pool, the nanoseq output is named by the barcode
                    else:
                        if int(barcode) < 10:
                            barcode = '0' + barcode  # Correction for guppy output filenames like barcode01.fastq.gz
                        anglerfish_content += sample_name + ',' + run_type + ',' + illumina_barcode + ',' + os.path.join(fastq_location, 'barcode' + barcode + '.fastq.gz') + '\n'
        with open(self.nanoseq_sample_sheet, 'w') as f:
            f.write(nanoseq_content)
        if anglerfish_content:
            with open(self.anglerfish_sample_sheet, 'w') as f:
                f.write(anglerfish_content)

    def check_exit_status(self, status_file):
        """Read pipeline exit status file and return True if 0, False if anything else"""
        with open(status_file, 'r') as f:
            exit_status = f.readline().strip()
        return exit_status == '0'

    def start_anglerfish(self):
        """Start Anglerfish."""
        os.makedirs(self.anglerfish_dir)
        anglerfish_command = ('anglerfish' #TODO: check that this is still correct
                            + ' --samplesheet ' + self.anglerfish_sample_sheet
                            + ' --out_fastq ' + self.anglerfish_dir
                            + ' --threads 2'
                            + ' --skip_demux; echo $? > .exitcode_for_anglerfish')
        try:
            p_handle = subprocess.Popen(anglerfish_command, stdout=subprocess.PIPE, shell=True, cwd=self.run_dir)
            logger.info('Started Anglerfish for run {} using: {}'.format(self.run_dir, anglerfish_command))
        except subprocess.CalledProcessError:
            logger.warn('An error occurred while starting the Anglerfish for run {}. '
                        'Please check the logfile for info.'.format(self.run_dir))

    def copy_results_for_lims(self):
        """Find results and copy to lims directory."""
        year_processed = self.run_id[0:4]
        lims_result_file = os.path.join(CONFIG.get('nanopore_analysis').get('minion_qc_run').get('lims_results_dir'),
                                        year_processed, 'anglerfish_stats_' + self.flowcell_id + '.txt')
        anglerfish_results = self._find_anglerfish_results()
        try:
            shutil.copyfile(anglerfish_results, lims_result_file)
            return True
        except TypeError as e:
            logger.warn('An error occurred while copying the Anglerfish results for {} to lims: {}'.format(self.run_id, e))
            return False

    def _find_anglerfish_results(self):
        """Return location of Anglerfish results."""
        results_file = ''
        for sub_dir in os.listdir(self.anglerfish_dir):
            if 'anglerfish_stats.txt' in os.listdir(os.path.join(self.anglerfish_dir, sub_dir)):
                results_file = os.path.join(self.anglerfish_dir, sub_dir, 'anglerfish_stats.txt')
                return results_file
        if not results_file:
            logger.warn('Could not find any Anglerfish results in {}'.format(self.anglerfish_dir))

class MinIONdelivery(Nanopore):
    """Minion delivery run"""
    def __init__(self, run_dir):
        super(MinIONdelivery, self).__init__(run_dir)
        self.transfer_details = CONFIG.get('nanopore_analysis').get('minion_delivery_run').get('transfer')
        self.transfer_log = self.transfer_details.get('transfer_file')
        self.archive_dir = self.transfer_details.get('finished_dir')
    
    def dump_path(self):
        """Dump path to run to a file that can be
        used when uploading stats to statusdb from preproc."""
        new_file = os.path.join(self.run_dir, 'run_path.txt')
        proj, sample, run = self.run_dir.split('/')[-3:]
        path_to_write = os.path.join(proj, sample, run)
        with open(new_file, 'w') as f:
            f.write(path_to_write)
    
    def write_finished_indicator(self):
        """Write a hidden file to indicate 
        when the finial rsync is finished."""
        new_file = os.path.join(self.run_dir, '.sync_finished')
        pathlib.Path(new_file).touch()
        return new_file