
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
	disk_space = status.get_disk_space()
	for key, value in disk_space.iteritems():
		print key, '	:	', value
	status.update_google_docs(disk_space)
