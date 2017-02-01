
""" CLI for the bioinfo subcommand
"""
import click
import taca.utils.bioinfo_tab as bt

@click.group()
def bioinfo_deliveries():
	""" Update statusdb with information about FC entry point """
	pass

# bioinfo subcommands
@bioinfo_deliveries.command()
@click.argument('rundir')
def updaterun(rundir):
    """saves the bioinfo data to statusdb"""
    bt.update_statusdb(rundir)

@bioinfo_deliveries.command()
def update():
    """saves the bioinfo data of everything that can be found to statusdb"""
    bt.collect_runs()

@bioinfo_deliveries.command()
@click.argument('runid')
@click.option('-p','--project', is_flag=False, help='Fail run for the specified project')
def fail_run(runid, project=None):
    """updates the status of the specified run to 'Failed'.
    Example of RUNID: 170113_ST-E00269_0163_BHCVH7ALXX"""
    bt.fail_run(runid, project)