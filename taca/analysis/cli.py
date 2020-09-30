"""CLI for the analysis subcommand."""
import click
from taca.analysis import analysis
from taca.analysis import analysis_nanopore

@click.group()
def analysis():
    """Analysis methods entry point."""
    pass

# Illumina analysis subcommands
@analysis.command()
@click.option('-r', '--run', type=click.Path(exists=True), default=None,
				 help='Demultiplex only a particular run')
@click.option('--force', is_flag=True, help='If specified, always tranfers the runs, even if they fail QC. Mail is sent anyway' )

def demultiplex(run, force):
	"""Demultiplex and transfer all runs present in the data directories."""
	analysis.run_preprocessing(run, force_trasfer=force)

@analysis.command()
@click.option('--runfolder-project', is_flag=False, help='Project ID for runfolder transfer')
@click.option('--exclude-lane', default='', help='Lanes to exclude separated by comma')
@click.argument('rundir')

def transfer(rundir, runfolder_project, exclude_lane):
    """Transfers the run without qc."""
    if not runfolder_project:
        analysis.transfer_run(rundir)
    else:
        analysis.transfer_runfolder(rundir, pid=runfolder_project, exclude_lane=exclude_lane)

@analysis.command()
@click.argument('rundir')
def updatedb(rundir):
    """Save the run to statusdb."""
    analysis.upload_to_statusdb(rundir)

# Nanopore analysis subcommans
@analysis.command()
@click.option('-r', '--run', type=click.Path(exists=True), default=None,
              help='Demultiplex only a particular run')
@click.option('--nanoseq_sample_sheet', type=click.Path(exists=True), default=None,
              help='Sample sheet for running nanoseq')
@click.option('--anglerfish_sample_sheet', type=click.Path(exists=True), default=None,
              help='Sample sheet for running anglerfish. Also requires --nanoseq_sample_sheet')

def demultiplex_nanopore(run, nanoseq_sample_sheet, anglerfish_sample_sheet):
    """Analyse and transfer all runs present in the data directories.
    Assumes QC run per default. Use --nanoseq_sample_sheet without --anglerfish_sample_sheet
    to manually start non-QC runs."""
    if anglerfish_sample_sheet and not nanoseq_sample_sheet:
        print('ERROR: Please specify --nanoseq_sample_sheet when using --anglerfish_sample_sheet')
        return
    analysis_nanopore.run_preprocessing(run, nanoseq_sample_sheet, anglerfish_sample_sheet)
