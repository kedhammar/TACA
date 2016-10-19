""" CLI for the storage subcommand
"""
import click
from taca.cleanup import cleanup as cln
from taca.utils import misc

@click.group()
@click.pass_context
def cleanup(ctx):
	""" Cleaning up servers - management methods and utilities """
	pass

# cleanup subcommands
@cleanup.command()
@click.option('-d', '--days', type=click.IntRange(min=1),
              help="Days to consider as thershold, should not be combined with option '--hours'")
@click.option('-h', '--hours', type=click.IntRange(min=1),
              help="Hours to consider as thershold, should not be combined with option '--days'")
@click.pass_context
def nas(ctx, days, hours):
    """ Do appropriate cleanup on NAS"""
    seconds = misc.to_seconds(days, hours)
    cln.cleanup_nas(seconds)

@cleanup.command()
@click.option('-d', '--days', type=click.IntRange(min=1),
              help="Days to consider as thershold, should not be combined with option '--hours'")
@click.option('-h', '--hours', type=click.IntRange(min=1),
              help="Hours to consider as thershold, should not be combined with option '--days'")
@click.pass_context
def preproc(ctx, days, hours):
    """ Do appropriate cleanup on PREPROC"""
    seconds = misc.to_seconds(days, hours)
    cln.cleanup_processing(seconds)

@cleanup.command()
@click.option('--site', type=click.Choice(['illumina','analysis','archive']),
              required=True, help="Which data to be cleaned on 'milou' server")
@click.option('-d', '--days', type=click.IntRange(min=1), help="Days to consider as thershold")
@click.option('-n','--dry-run', is_flag=True, help='Perform dry run i.e. Executes nothing but log')
@click.pass_context
def milou(ctx, site, days, dry_run):
    """ Do appropriate cleanup on MILOU"""
    seconds = misc.to_seconds(days)
    cln.cleanup_milou(site, seconds, dry_run)

@cleanup.command()
@click.option('--days_fastq', type=click.IntRange(min=1),
              help="Days to consider as thershold for removing 'fastq' files")
@click.option('--days_analysis', type=click.IntRange(min=1),
              help="Days to consider as thershold for removing analysis data")
@click.option('--only_fastq', is_flag=True, help="Clean only fastq data in 'irma'")
@click.option('--only_analysis', is_flag=True, help="Clean only analysis data in 'irma'")
@click.option('-n','--dry_run', is_flag=True, help='Perform dry run i.e. Executes nothing but log')
@click.pass_context
def irma(ctx, days_fastq, days_analysis, only_fastq, only_analysis, dry_run):
    """ Do appropriate cleanup on IRMA"""
    config_file = ctx.parent.parent.params['config_file'].name
    if only_fastq and only_analysis:
        raise SystemExit("ERROR: Both option 'only_fastq' and 'only_analysis' is given, should only give either one")
    if not days_fastq and not only_analysis:
        raise SystemExit("ERROR: 'days_fastq' is not given while not selecting 'only_analysis' option")
    if not days_analysis and not only_fastq:
        raise SystemExit("ERROR: 'days_analysis' is not given while not selecting 'only_fastq' option")
    cln.cleanup_irma(config_file, days_fastq, days_analysis, only_fastq, only_analysis, dry_run)
