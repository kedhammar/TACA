
""" CLI for the testing commands
"""
import click
import taca.testing.create_uppmax_like_env as createupp


@click.group()
def uppmax_env():
    """ Create a local set of folders that resembles the uppmax-ngi env. Creates config file for ngi_pipeline, taca, and taca ngi-pipeline. Only a minimal taca config is needed (statusdb and log) """
    pass

@uppmax_env.command()
@click.option('-p', '--projects', type=int, default=30, help='number of projects to be extracted from statusdb')

def create(projects):
    """creates a uppmax like env 
    """
    createupp.create(projects)
