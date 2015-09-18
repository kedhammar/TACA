
import click
import logging
import os

from taca.server_status import server_status as status
from taca.utils.config import CONFIG

@click.group()
def server_status():
	""" Monitor server status """
	if not CONFIG.get('server_status', ''):
		raise RuntimeError("Configuration missing required entries: server_status")



# server status subcommands
@server_status.command()
@click.option('--credentials', type=click.Path(exists=True), default=os.path.join(os.path.dirname(__file__), 'gdocs_credentials.json'),
				 help='Path to google credentials file')
def disk_space(credentials):
	""" Checks the available space on all the disks
	"""
	disk_space = status.get_disk_space()
	status.update_google_docs(disk_space, credentials)
