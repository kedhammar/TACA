"""CLI for the analysis subcommand."""
import click

from taca.analysis import analysis as an
from taca.analysis import analysis_nanopore

@click.group()
def analysis():
    """Analysis methods entry point."""
    pass

# Illumina analysis subcommands
@analysis.command()
@click.option('-r', '--run', type=click.Path(exists=True), default=None,
				 help='Demultiplex only a particular run')

def demultiplex(run):
	"""Demultiplex and transfer all runs present in the data directories."""
	an.run_preprocessing(run)

@analysis.command()
@click.option('--runfolder-project', is_flag=False, help='Project ID for runfolder transfer')
@click.option('--exclude-lane', default='', help='Lanes to exclude separated by comma')
@click.argument('rundir')

def transfer(rundir, runfolder_project, exclude_lane):
    """Transfers the run without qc."""
    if not runfolder_project:
        an.transfer_run(rundir)
    else:
        an.transfer_runfolder(rundir, pid=runfolder_project, exclude_lane=exclude_lane)

@analysis.command()
@click.argument('rundir')
def updatedb(rundir):
    """Save the run to statusdb."""
    an.upload_to_statusdb(rundir)

# Nanopore analysis subcommands
@analysis.command()
@click.option('--runtype', type=click.Choice(['qc', 'delivery']),
              help='Required. Specify if MinION run is for QC or delivery.')
@click.option('-r', '--run', type=click.Path(exists=True), default=None,
              help='Process only a particular run')
@click.option('--nanoseq_sample_sheet', type=click.Path(exists=True), default=None,
              help='Manually edited sample sheet for running nanoseq')
@click.option('--anglerfish_sample_sheet', type=click.Path(exists=True), default=None,
              help='Manually edited sample sheet for running anglerfish')

def minion(runtype, run, nanoseq_sample_sheet, anglerfish_sample_sheet):
    """Process MinION QC runs
    """
    if runtype == 'qc':
        analysis_nanopore.process_minion_qc_runs(run, nanoseq_sample_sheet, anglerfish_sample_sheet)
    elif runtype == 'delivery':
        analysis_nanopore.process_minion_delivery_runs(run)
    else:
        print('Please specify the MinION runtype (qc or delivery)')

@analysis.command()
@click.option('-r', '--run', type=click.Path(exists=True), default=None,
              help='Transfer only a particular run')

def ont_transfer(run):
    """Transfer runs present in the data directories to HPC cluster.
    """
    analysis_nanopore.transfer_finished(run)