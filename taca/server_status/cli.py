import click
import logging

from taca.server_status import server_status as status
from taca.utils.config import CONFIG
from taca.server_status import cronjobs as cj # to avoid similar names with command, otherwise exception


@click.group(name='server_status')
def server_status():
    """ Monitor server status """

# server status subcommands
@server_status.command()
@click.option('--statusdb', is_flag=True, help="Update the statusdb")
def nases(statusdb):
    """ Checks the available space on all the nases
    """
    if not CONFIG.get('server_status', ''):
        logging.warning("Configuration missing required entries: server_status")
    disk_space = status.get_nases_disk_space()
    if statusdb:
        status.update_status_db(disk_space, server_type='nas')

@server_status.command()
def cronjobs():
    """ Monitors cronjobs and updates statusdb
    """
    cj.update_cronjob_db()

@server_status.command()
def monitor_promethion():
    """ Checks the status of PromethION and if ngi-nas is mounted
    """
    if not CONFIG.get('promethion_status', ''):
        logging.warning("Configuration missing required entries: server_status")
    promethion_status = status.check_promethion_status()
    if promethion_status:
        logging.info("No issues encountered with the PromethION")
    else:
        logging.warning("An issue with the PromethION was encountered. Operator has been notified by email.")