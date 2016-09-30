
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
    createupp.create(projects, ngi_config)
