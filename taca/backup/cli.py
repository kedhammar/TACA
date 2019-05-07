""" CLI for the backup subcommand
"""
import click
from taca.backup.backup import backup_utils as bkut

@click.group()
@click.pass_context
def backup(ctx):
	""" Backup management methods and utilities """
	pass


@backup.command()
@click.option('-r', '--run', type=click.Path(exists=True), help="A run (directory or a zipped archive) to be encrypted")
@click.option('-f', '--force', is_flag=True, help="Ignore the checks and just try encryption. USE IT WITH CAUTION.")
@click.pass_context
def encrypt(ctx, run, force):
    bkut.encrypt_runs(run, force)

@backup.command(name='put_data')
@click.option('-r', '--run', type=click.Path(exists=True), help="A run name (without extension) to be sent to PDC")
@click.pass_context
def put_data(ctx, run):
    bkut.pdc_put(run)

@backup.command(name='get_data')
@click.option('-r', '--run', required=True, help="A run name (without extension) to download from PDC")
@click.option('-o', '--outdir', type=click.Path(exists=True, file_okay=False, writable=True),
              help="Optional directory name to save the downloaded file. Directory should exist")
@click.pass_context
def get_data(ctx, run, outdir):
    ## W I P ##
    raise NotImplementedError

@backup.command()
@click.option('-r', '--run', required=True, type=click.Path(exists=True, dir_okay=False), help="A encripted run file")
@click.option('-k', '--key', required=True, help="Key file to be used for decryption")
@click.option('-p', '--password', help="To pass decryption passphrase via command line")
@click.pass_context
def decrypt(ctx, run, key, password):
    ## W I P ##
    raise NotImplementedError
