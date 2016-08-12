""" CLI for the backup subcommand
"""
import click
from taca.backup import backup as bk
from taca.backup.backup_util import backup_utils as bkut
        
@click.group()
@click.pass_context
def backup(ctx):
	""" Backup management methods and utilities """
	pass


@backup.command()
@click.option('-r', '--run', type=click.Path(exists=True), help="A run (directory or a zipped archive) to be encrypted")
@click.option('-f', '--force', is_flag=True, help="Ignore the checks and just try encryption. USE IT WITH CAUTION.")
@click.option('-m', '--method_call', is_flag=True)
@click.option('-c', '--class_call', is_flag=True)
@click.pass_context
def encrypt(ctx, run, force, method_call, class_call):
    if method_call:
        bk.encrypt_data(run, force)
    if class_call
        bkut.encrypt_runs(run, force)

@backup.command()
@click.option('-r', '--run', type=click.Path(exists=True), help="A run name (without extension) to be sent to PDC")
@click.pass_context
def put_data(ctx, run):
    ## W I P ##
    raise NotImplementedError

@backup.command()
@click.option('-r', '--run', required=True, help="A run name (without extension) to download from PDC")
@click.pass_context
def get_data(ctx, run):
    ## W I P ##
    raise NotImplementedError

@backup.command()
@click.option('-r', '--run', required=True, help="A run name (without extension) to download from PDC")
@click.pass_context
def decrypt(ctx, run):
    ## W I P ##
    raise NotImplementedError
