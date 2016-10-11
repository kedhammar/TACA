
""" CLI for the testing commands
"""
import os
import click
import taca.testing.create_uppmax_like_env as createupp

@click.group()
def uppmax_env():
    """ Create a local set of folders that resembles the uppmax-ngi env. Creates config file for ngi_pipeline, taca, and taca ngi-pipeline. Only a minimal taca config is needed (statusdb and log) 
        The condig file (in general saved in variable NGI_CONFIG needs to looks something similar to:
        
        \b
        environment:
            project_id: ngi1234 #CAN BE ANYTHING
            ngi_scripts_dir: /Users/vezzi/opt/ngi_pipeline/scripts #CAN BE ANYTHING
            conda_env: TACA #CAN BE ANYTHING
        flowcell_inbox:
            - /Users/vezzi/opt/uppmax_env/incoming/ #NEEDS TO EXISTS
        analysis:
            best_practice_analysis:
                whole_genome_reseq:
                    analysis_engine: ngi_pipeline.engines.piper_ngi
                IGN:
                    analysis_engine: ngi_pipeline.engines.piper_ngi
                    
                qc:
                
                    analysis_engine: ngi_pipeline.engines.qc_ngi
                    
            base_root: /Users/vezzi/opt/  #NEEDS TO EXISTS
            sthlm_root: uppmax_env        #NEEDS TO EXISTS
            top_dir: nobackup/NGI         #NEEDS TO EXISTS
            upps_root: nothing            #CAN BE ANYTHING
        logging:
            log_file: "/Users/vezzi/opt/log/ngi_pipeline.log" #NEEDS TO BE REAL
  
        \b
        The requested project will be divided into the following sets:
          - 2/3 will be selected among the projects with application equeal to 'WG re-seq'. These will be divided up in:
            - 1/4: closed more than 3 months ago
            - 1/4: closed more than 1 month ago, less than 3 months
            - 1/4: closed less than 1 month ago
            - 1/4: open
          - 1/3 will be selected amonf the projects with application different from 'WG re-seq':
            - 1/4: closed more than 3 months ago
            - 1/4: closed more than 1 month ago, less than 3 months
            - 1/4: closed less than 1 month ago
            - 1/4: open
   
     """
    pass

@uppmax_env.command()
@click.option('-p', '--projects', type=int, default=30, help='number of projects to be extracted from statusdb')
@click.option('-nc', '--ngi-config', type=str,  default=os.environ.get('NGI_CONFIG') , help='path to ngi configuration file (expected in variable NGI_CONFIG)')


def create(projects, ngi_config):
    """creates a uppmax like env 
    """
    if which("ngi_pipeline_start.py"):
        createupp.create(projects, ngi_config)
    else:
        print "ERROR: ngi_pipeline_start.py needs to be available and properly installed"


def which(file):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
                return True
    return False