""" CLI for the analysis subcommand
"""
import click
from taca.analysis import analysis
from taca.analysis import analysis_nanopore

@click.group()
def analysis():
	""" Analysis methods entry point """
	pass

# Illumina analysis subcommands
@analysis.command()
@click.option('-r', '--run', type=click.Path(exists=True), default=None,
				 help='Demultiplex only a particular run')
@click.option('--force', is_flag=True, help='If specified, always tranfers the runs, even if they fail QC. Mail is sent anyway' )

def demultiplex(run, force):
	""" Demultiplex and transfer all runs present in the data directories
	"""
	analysis.run_preprocessing(run, force_trasfer=force)

@analysis.command()
@click.option('-a','--analysis', is_flag=False, help='Trigger the analysis for the transferred flowcell')
@click.option('--runfolder-project', is_flag=False, help='Project ID for runfolder transfer')
@click.argument('rundir')

def transfer(rundir, analysis, runfolder_project):
    """ Transfer the run without qc"""
    if not runfolder_project:
        analysis.transfer_run(rundir, analysis=analysis)
    else:
        analysis.transfer_runfolder(rundir, pid=runfolder_project)

@analysis.command()
@click.argument('rundir')
def updatedb(rundir):
    """ Save the run to statusdb"""
    analysis.upload_to_statusdb(rundir)

# Nanopore analysis subcommans
@analysis.command()
@click.option('-r', '--run', type=click.Path(exists=True), default=None,
              help='Demultiplex only a particular run')

def demultiplex_nanopore(run):
    """ Basecall, demultiplex and transfer all runs present in the data directories"""
    analysis_nanopore.run_preprocessing(run)
