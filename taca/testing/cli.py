
""" CLI for the testing commands
"""
import os
import click
import taca.testing.create_uppmax_like_env as createupp

@click.group()
def uppmax_env():
    """ Create a local set of folders that resembles the uppmax-ngi env. Creates config file for ngi_pipeline, taca, and taca ngi-pipeline. Only a minimal taca config is needed (statusdb and log) """
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
        sys.exit("ERROR: ngi_pipeline_start.py needs to be available and properly installed")


def which(file):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, file)):
                return True
    return False