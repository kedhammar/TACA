""" CLI for the storage subcommand
"""
import click
from taca.storage import storage as st
from taca.utils import misc
        
@click.group()
@click.pass_context
def storage(ctx):
	""" Storage management methods and utilities """
	pass

# Storage subcommands
@storage.command()
@click.option('-d', '--days', type=click.IntRange(min=1),
              help="Days to consider as thershold, should not be combined with option '--hours'")
@click.option('-h', '--hours', type=click.IntRange(min=1),
              help="Hours to consider as thershold, should not be combined with option '--days'")
@click.option('-s','--site', type=click.Choice(['archive','illumina','analysis','nas','processing-server']),
              required=True, help='Site to perform cleanup')
@click.option('-n','--dry-run', is_flag=True, help='Perform dry run i.e. Executes nothing but log')
@click.pass_context
def cleanup(ctx, days, hours, site, dry_run):
    """ Do appropriate cleanup on the given site i.e. NAS/processing servers/UPPMAX """
    seconds = misc.to_seconds(days, hours)
    if site == 'nas':
        st.cleanup_nas(seconds)
    if site == 'processing-server':
        st.cleanup_processing(seconds)
    if site in ['illumina','analysis','archive']:
        st.cleanup_uppmax(site, seconds, dry_run)
