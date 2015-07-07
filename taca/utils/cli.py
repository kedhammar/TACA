
""" CLI for the bioinfo subcommand
"""
import click
import taca.utils.bioinfo_tab as bt


@click.group()
def bioinfo():
	""" Analysis methods entry point """
	pass

# bioinfo subcommands
@bioinfo.command()
@click.argument('rundir')
def updatedb(rundir):
    """saves the bioinfo data to statusdb"""
    bt.update_statusdb(rundir)

@bioinfo.command()
def update():
    """saves the bioinfo data of everything that can be found to statusdb"""
    bt.collect_runs()
