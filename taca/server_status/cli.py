import click
import logging
import os

from taca.server_status import server_status as status
from taca.utils.config import CONFIG
from taca.server_status import cronjobs as cj # to avoid similar names with command, otherwise exception


@click.group(name='server_status')
def server_status():
    """ Monitor server status """
    if not CONFIG.get('server_status', ''):
        logging.warning("Configuration missing required entries: server_status")

# server status subcommands
@server_status.command()
@click.option('--statusdb', is_flag=True, help="Update the statusdb")
def nases(statusdb):
    """ Checks the available space on all the nases
    """
    disk_space = status.get_nases_disk_space()
    if statusdb:
        status.update_status_db(disk_space, server_type='nas')

@server_status.command()
def cronjobs():
    """ Monitors cronjobs and updates statusdb
    """
    cj.update_cronjob_db()
