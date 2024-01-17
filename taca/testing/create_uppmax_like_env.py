""" Load and parse configuration file."""

import datetime
import logging
import os
import random
import subprocess
from io import open

from dateutil.relativedelta import relativedelta

from taca.utils import config as conf
from taca.utils import filesystem as fs
from taca.utils import statusdb
from taca.utils.config import CONFIG

logger = logging.getLogger(__name__)


def create_version_report(path):
    # Creates the file version_report.txt for stuff run ngi_pipeline
    with open(os.path.join(path, 'version_report.txt'), 'w') as VERSION_REPORT:
        VERSION_REPORT.write('******\n')
        VERSION_REPORT.write('README\n')
        VERSION_REPORT.write('******\n')
        VERSION_REPORT.write('\n')
        VERSION_REPORT.write('Data has been aligned to to the reference using bwa. The raw alignments have then been deduplicated, recalibrated and cleaned using GATK. Quality control information was gathered using Qualimap. SNVs and indels have been called using the HaplotypeCaller. These variants were then funcionally annotated using snpEff. The pipeline used was Piper, see below for more information.\n')
        VERSION_REPORT.write('\n')
        VERSION_REPORT.write('The versions of programs and references used:\n')
        VERSION_REPORT.write('piper: unknown\n')
        VERSION_REPORT.write('bwa: 0.7.12\n')
        VERSION_REPORT.write('samtools: 0.1.19\n')
        VERSION_REPORT.write('qualimap: v2.2\n')
        VERSION_REPORT.write('snpEff: 4.1\n')
        VERSION_REPORT.write('snpEff reference: GRCh37.75\n')
        VERSION_REPORT.write('gatk: 3.3-0-geee94ec\n')
        VERSION_REPORT.write('\n')
        VERSION_REPORT.write('reference: human_g1k_v37.fasta\n')
        VERSION_REPORT.write('db_snp: gatk-bundle/2.8\n')
        VERSION_REPORT.write('hapmap: gatk-bundle/2.8\n')
        VERSION_REPORT.write('omni: gatk-bundle/2.8\n')
        VERSION_REPORT.write('1000G_indels: gatk-bundle/2.8\n')
        VERSION_REPORT.write('Mills_and_1000G_golden_standard_indels: gatk-bundle/2.8\n')
        VERSION_REPORT.write('\n')
        VERSION_REPORT.write('indel resource file: {Mills_and_1000G_gold_standard.indels.b37.vcf version: gatk-bundle/2.8}\n')
        VERSION_REPORT.write('indel resource file: {1000G_phase1.indels.b37.vcf version: gatk-bundle/2.8}\n')
        VERSION_REPORT.write('\n')
        VERSION_REPORT.write('piper\n')
        VERSION_REPORT.write('-----\n')
        VERSION_REPORT.write('Piper is a pipeline system developed and maintained at the National Genomics Infrastructure build on top of GATK Queue. For more information and the source code visit: www.github.com/NationalGenomicsInfrastructure/piper\n')

def create_FC(incoming_dir, run_name, samplesheet, fastq_1 = None, fastq_2=None ):
    # Create something like 160217_ST-E00201_0063_AHJHNYCCXX
    path_to_fc = os.path.join(incoming_dir, run_name)
    if os.path.exists(path_to_fc):
        # This FC exists, skip it
        return
    fs.create_folder(path_to_fc)
    fs.touch(os.path.join(path_to_fc, 'RTAComplete.txt'))
    # Create folder Demultiplexing
    fs.create_folder(os.path.join(path_to_fc, 'Demultiplexing'))
    # Create folder Demultiplexing/Reports
    fs.create_folder(os.path.join(path_to_fc, 'Demultiplexing', 'Reports'))
    # Create folder Demultiplexing/Stats
    fs.create_folder(os.path.join(path_to_fc, 'Demultiplexing', 'Stats'))
    # Memorise SampleSheet stats
    header = []
    for key in samplesheet[0]:
        header.append(key)
    counter = 1
    current_lane = ''
    for line in samplesheet:
        project_name = line.get('Sample_Project', line.get('Project', ''))
        lane = line['Lane']
        if current_lane == '':
            current_lane = lane
        elif current_lane != lane:
            counter = 1
            current_lane = lane
        sample_id = line.get('SampleID', line.get('Sample_ID', ''))
        sample_name = line.get('SampleName', line.get('Sample_Name', ''))
        # Create dir structure
        fs.create_folder(os.path.join(path_to_fc, 'Demultiplexing', project_name, sample_id))
        # Now create the data
        fastq_1_dest = f'{sample_name}_S{counter}_L00{lane}_R1_001.fastq.gz'
        fastq_2_dest = f'{sample_name}_S{counter}_L00{lane}_R2_001.fastq.gz'
        counter += 1
        if fastq_1 is None:
            fs.touch(os.path.join(path_to_fc, 'Demultiplexing', project_name,
                                  sample_id, fastq_1_dest))
            fs.touch(os.path.join(path_to_fc, 'Demultiplexing', project_name,
                                  sample_id, fastq_2_dest))
        else:
            fs.do_symlink(fastq_1, os.path.join(path_to_fc, 'Demultiplexing',
                                                project_name, sample_id, fastq_1_dest))
            fs.do_symlink(fastq_2, os.path.join(path_to_fc, 'Demultiplexing',
                                                project_name, sample_id, fastq_2_dest))

    with open(os.path.join(path_to_fc, 'SampleSheet.csv'), 'w') as Samplesheet_file:
        Samplesheet_file.write('[Header]\n')
        Samplesheet_file.write('Date,2016-03-29\n')
        Samplesheet_file.write('Investigator Name,Christian Natanaelsson\n')
        Samplesheet_file.write('[Data]\n')
        for key in header:
             Samplesheet_file.write(f'{key},')
        Samplesheet_file.write('\n')
        for line in samplesheet:
            for key in header:
                Samplesheet_file.write(f'{line[key]},')
            Samplesheet_file.write('\n')

def create_uppmax_env(ngi_config):
    paths = {}
    if 'analysis' not in ngi_config:
        sys.exit('ERROR: analysis must be a field of NGI_CONFIG.')
    try:
        base_root = ngi_config['analysis']['base_root']
        paths['base_root'] = base_root
        sthlm_root = ngi_config['analysis']['sthlm_root']
        paths['sthlm_root'] = sthlm_root
        top_dir = ngi_config['analysis']['top_dir']
        paths['top_dir'] = top_dir
    except KeyError as e:
        raise SystemExit(f'Config file is missing the key {str(e)}, make sure it have all required information')
    if 'environment' not in ngi_config:
        sys.exit('ERROR: environment must be a field of NGI_CONFIG.')
    try:
        # Get base root
        flowcell_inboxes = ngi_config['environment']['flowcell_inbox']
        flowcell_inbox = flowcell_inboxes[0] # I assume there is only one
        paths['flowcell_inbox'] = flowcell_inbox
    except ValueError as e:
        sys.exit(f'key error, flowcell_inbox not found in "{ngi_config}": {e}')
    # Now I need to create the folders for this
    if not os.path.exists(base_root):
        sys.exit(f'base_root needs to exists: {base_root}')
    fs.create_folder(flowcell_inbox)
    if sthlm_root is None:
        path_to_analysis = os.path.join(base_root, top_dir)
    else:
        path_to_analysis = os.path.join(base_root, sthlm_root, top_dir)
    fs.create_folder(path_to_analysis)
    return paths

def produce_analysis_qc_ngi(ngi_config, project_id):
    analysis_dir = os.path.join(ngi_config['analysis']['base_root'],
                                ngi_config['analysis']['sthlm_root'],
                                ngi_config['analysis']['top_dir'],
                                'ANALYSIS', project_id)
    data_dir = os.path.join(ngi_config['analysis']['base_root'],
                            ngi_config['analysis']['sthlm_root'],
                            ngi_config['analysis']['top_dir'],
                            'DATA', project_id)

    qc_ngi_dir = os.path.join(analysis_dir, 'qc_ngi')
    fs.create_folder(qc_ngi_dir)
    for sample_id in os.listdir(data_dir):
        sample_dir_qc = os.path.join(qc_ngi_dir, sample_id)
        fs.create_folder(sample_dir_qc)
        fastqc_dir = os.path.join(sample_dir_qc, 'fastqc')
        fs.create_folder(fastqc_dir)
        fastq_screen_dir  = os.path.join(sample_dir_qc, 'fastq_screen')
        fs.create_folder(fastq_screen_dir)
        # Do not create more than this...

def produce_analysis_piper(ngi_config, project_id):
    # Create piper_ngi
    analysis_dir = os.path.join(ngi_config['analysis']['base_root'],
                                ngi_config['analysis']['sthlm_root'],
                                ngi_config['analysis']['top_dir'],
                                'ANALYSIS', project_id)
    data_dir = os.path.join(ngi_config['analysis']['base_root'],
                            ngi_config['analysis']['sthlm_root'],
                            ngi_config['analysis']['top_dir'],
                            'DATA', project_id)

    piper_ngi_dir = os.path.join(analysis_dir, 'piper_ngi')
    fs.create_folder(piper_ngi_dir)
    piper_dirs = ['01_raw_alignments',
                  '02_preliminary_alignment_qc',
                  '03_genotype_concordance',
                  '04_merged_aligments',
                  '05_processed_alignments',
                  '06_final_alignment_qc',
                  '07_variant_calls',
                  '08_misc']
    for piper_dir in piper_dirs:
        current_dir =  os.path.join(piper_ngi_dir, piper_dir)
        fs.create_folder(current_dir)
        if piper_dir == '05_processed_alignments':
            for sample_id in os.listdir(data_dir):
                bam_file = f'{sample_id}.clean.dedup.bam'
                fs.touch(os.path.join(current_dir, bam_file))
        if piper_dir == '07_variant_calls':
            for sample_id in os.listdir(data_dir):
                vcf_file = f'{sample_id}.clean.dedup.recal.bam.raw.indel.vcf.gz'
                fs.touch(os.path.join(current_dir, vcf_file))
    current_dir = os.path.join(piper_ngi_dir, 'sbatch')
    fs.create_folder(current_dir)
    current_dir = os.path.join(piper_ngi_dir, 'setup_xml_files')
    fs.create_folder(current_dir)
    current_dir = os.path.join(piper_ngi_dir, 'logs')
    fs.create_folder(current_dir)
    create_version_report(current_dir)

def select_random_projects(projects_in, num_proj, application, projects_out, label):
    chosen_projects = 0
    iterations = 0 # Safe guard to avoid infinite loops
    application_not_in_other = ['WG re-seq']
    while chosen_projects != num_proj and iterations < 4*len(projects_in):
        iterations += 1
        selected_proj = random.choice(list(projects_in.keys()))
        # Check if I have already picked up this element
        already_chosen = False
        for project_pair in projects_out:
            if selected_proj == project_pair[0]:
                already_chosen = True
        if already_chosen:
            continue # I am reprocessing an element I already saw. I skip it. iterations will avoid infinite loops
        proj_value = projects_in[selected_proj]
        if application == 'other':
            # In this case everything expcept
            if proj_value['application'] not in application_not_in_other:
                # I select this one
                projects_out.append([selected_proj, label])
                chosen_projects += 1
        elif application == proj_value['application']:
            # I select this one
            projects_out.append([selected_proj, label])
            chosen_projects += 1

def create(projects, ngi_config_file, fastq_1, fastq_2):
    statusdb_conf = CONFIG.get('statusdb')
    if statusdb_conf is None:
        logger.error('No statusdb field in taca configuration file')
        return 1
    if 'dev' not in statusdb_conf['url']:
        logger.error('url for status db is {}, but dev must be specified in this case'.format(statusdb_conf['url']))
    couch_connection = statusdb.StatusdbSession(statusdb_conf).connection
    projectsDB = couch_connection['projects']
    project_summary = projectsDB.view('project/summary')
    projects_closed_more_than_three_months = {}
    projects_closed_more_than_one_month_less_than_three = {}
    projects_closed_less_than_one_month = {}
    projects_opened = {}
    current_date = datetime.datetime.today()
    date_limit_one_year = current_date - relativedelta(months=6) #yes yes I know.. but in this way i am sure all data in in xflocell_db
    date_limit_one_month = current_date - relativedelta(months=1)
    date_limit_three_month = current_date - relativedelta(months=3)
    for row in project_summary:
        project_id = row['key'][1]
        project_status = row['key'][0]
        if 'application' not in row['value']:
            continue
        if row['value']['no_samples'] > 50:
            continue # Skip large projects
        application = row['value']['application']
        if project_status == 'closed':
            if 'close_date' in row['value']:
                close_date = datetime.datetime.strptime(row['value']['close_date'], '%Y-%m-%d')
                if close_date > date_limit_one_year: # If the project has been closed after the date limit
                    if close_date >= date_limit_one_month:
                        projects_closed_less_than_one_month[project_id] = {'project_name': row['value']['project_name'],
                                                                           'application': application,
                                                                           'no_samples': row['value']['no_samples']}
                    elif close_date < date_limit_one_month and close_date >= date_limit_three_month:
                        projects_closed_more_than_one_month_less_than_three[project_id] = {'project_name': row['value']['project_name'],
                                                                                           'application': application,
                                                                                           'no_samples': row['value']['no_samples']}
                    elif close_date < date_limit_three_month:
                        projects_closed_more_than_three_months[project_id] = {'project_name': row['value']['project_name'],
                                                                              'application': application,
                                                                              'no_samples': row['value']['no_samples']}
        elif project_status == 'open':
            if 'lanes_sequenced' in row['value'] and row['value']['lanes_sequenced'] > 0:
                projects_opened[project_id] =  {'project_name': row['value']['project_name'],
                                                'application': application,
                                                'no_samples': row['value']['no_samples']}
        else:
            print(f'status {project_status}')
    ## Now I can parse the x_flowcell db to check what I can and cannot use
    whole_genome_projects = int(2*projects/3)
    projects_to_reproduce = []
    select_random_projects(projects_closed_more_than_three_months,
                           whole_genome_projects/4+1,
                           'WG re-seq',
                           projects_to_reproduce,
                           'WGreseq_tot_closed')
    select_random_projects(projects_closed_more_than_one_month_less_than_three,
                           whole_genome_projects/4+1,
                           'WG re-seq',
                           projects_to_reproduce,
                           'WGreseq_closed_clean_no_del')
    select_random_projects(projects_closed_less_than_one_month,
                           whole_genome_projects/4+1,
                           'WG re-seq',
                           projects_to_reproduce,
                           'WGreseq_closed_no_clean')
    select_random_projects(projects_opened,
                           whole_genome_projects/4+1,
                           'WG re-seq',
                           projects_to_reproduce,
                           'WGreseq_open')

    other_projects = int(projects/3)
    select_random_projects(projects_closed_more_than_three_months,
                           other_projects/4+1,
                           'other',
                           projects_to_reproduce,
                           'noWGreseq_tot_closed')
    select_random_projects(projects_closed_more_than_one_month_less_than_three,
                           other_projects/4+1,
                           'other',
                           projects_to_reproduce,
                           'noWGreseq_closed_clean_no_del')
    select_random_projects(projects_closed_less_than_one_month,
                           other_projects/4+1,
                           'other',
                           projects_to_reproduce,
                           'noWGreseq_closed_no_clean')
    select_random_projects(projects_opened,
                           other_projects/4+1,
                           'other',
                           projects_to_reproduce,
                           'noWGreseq_open')

    # Create ngi_pipeline enviroment
    print(f'#NGI_CONFIG varaible is {ngi_config_file}. This variable needs to be in the .bashrc file')
    print(f'NGI_CONFIG={ngi_config_file}')
    try:
        ngi_config = conf.load_config(ngi_config_file)
    except OSError as e:
        print(f'ERROR: {e.message}')
    # Create uppmax env
    paths = create_uppmax_env(ngi_config)

    print(f'#Going to reproduce {len(projects_to_reproduce)} projects (if this number is different from the one you specified.... trust me... do not worry')
    # Scan over x_flowcell and reproduce FCs
    flowcellDB = couch_connection['x_flowcells']
    reproduced_projects = {}
    for fc_doc in flowcellDB:
        try:
            samplesheet_csv = flowcellDB[fc_doc]['samplesheet_csv']
        except KeyError:
            continue # Parse only FC that have a samplesheet
        # Check if this FC contains one of the proejcts I need to replicate.
        projects_in_FC = set()
        if 'SampleName' in samplesheet_csv[0]:
            projects_in_FC = set([line['SampleName'].split('_')[0] for line in samplesheet_csv])
        else:
            projects_in_FC = set([line['Sample_Name'].split('_')[0] for line in samplesheet_csv])
        found = False
        for project_pair in projects_to_reproduce:
            project = project_pair[0]
            if project in projects_in_FC:
                # This FC needs to be created
                if not found:
                    # Create the FC only the first time I see a project belonging to it
                    create_FC(paths['flowcell_inbox'] , flowcellDB[fc_doc]['RunInfo']['Id'], samplesheet_csv, fastq_1, fastq_2)
                    found = True
                # But I keep track of all projects-run I need to organise
                if project not in reproduced_projects:
                    reproduced_projects[project] = []
                reproduced_projects[project].append(flowcellDB[fc_doc]['RunInfo']['Id'])
    print(f'#Reproduced {len(reproduced_projects)} project (if the numbers diffear do not worry, most likely we selected projects without runs)')
    for project in projects_to_reproduce:
        if project[0] in reproduced_projects:
            print(f'#  {project[0]}: {project[1]}')
    # Need to output the command to organise
    to_be_deleted = []
    for project in reproduced_projects:
        for FC in reproduced_projects[project]:
            print(f'Running: ngi_pipeline_start.py organize flowcell {FC} -p {project}')
            with open('ngi_pipeline_local.logs', 'w') as NGILOGS:
                return_value = subprocess.call(['ngi_pipeline_start.py',
                                                'organize',
                                                'flowcell',
                                                f'{FC}',
                                                '-p',
                                                f'{project}'],
                                               stdout=NGILOGS, stderr=NGILOGS)
            if return_value > 0:
                print(f'#project {project} not organised: have a look to the logs, but most likely this projec is not in charon')
                if project not in to_be_deleted:
                    to_be_deleted.append(project)

    for project in to_be_deleted:
        del reproduced_projects[project]

    # Create ANALYSIS --
    for project in projects_to_reproduce:
        if project[0] in reproduced_projects: # Only for projects that I know I have organised
            produce_analysis_qc_ngi(ngi_config, project[0])
            if project[1].startswith('WGreseq'):
                produce_analysis_piper(ngi_config, project[0])

    # Store in a file the results
    with open('projects.txt', 'w') as PROJECTS:
        for project in projects_to_reproduce:
            if project[0] in reproduced_projects:
                PROJECTS.write(f'{project[0]}:{project[1]}\n')
