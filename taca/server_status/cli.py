
import click
from taca.server_status import server_status as status

@click.group()
def server_status():
	""" Monitor server status """
	pass


# server status subcommands
@server_status.command()
def disk_space():
	""" Checks the available space on all the disks
	"""
	disk_space = status.disk_space()
	print disk_space
	return disk_space