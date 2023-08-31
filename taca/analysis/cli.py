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
@click.option('-r', '--run', type=click.Path(exists=True), default=None,
              help='Transfer only a particular run')
def ont_transfer(run):
    """Transfer runs present in the data directories to HPC cluster.
    """
    analysis_nanopore.transfer_finished(run)

@analysis.command()
@click.argument("run")
def ont_updatedb(run):
    """Update the database, regardless of run status"""
    analysis_nanopore.ont_updatedb_from_cli(run)