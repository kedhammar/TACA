
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
@click.option('--statusdb', is_flag=True, help="Update the statusdb")
def nases(statusdb):
	""" Checks the available space on all the nases
	"""
	disk_space = status.get_nases_disk_space()
	if statusdb:
		status.update_status_db(disk_space, server_type='nas')

#  must be run on uppmax, as no passwordless ssh to uppmax servers
@server_status.command()
@click.option('--disk-quota', is_flag=True, help="Check the available space on the disks")
@click.option('--cpu-hours', is_flag=True, help="Check the usage of CPU hours")
def uppmax(disk_quota, cpu_hours):
	"""
	Checks the quotas and cpu hours on the uppmax servers
	"""
	merged_results = {}
	if disk_quota:
		disk_quota_data = status.get_uppmax_quotas()
		merged_results.update(disk_quota_data)
	if cpu_hours:
		cpu_hours_data = status.get_uppmax_cpu_hours()
		if not merged_results:
			merged_results = cpu_hours_data
		else:
			for key in cpu_hours_data.keys():
				if key not in merged_results:
					merged_results[key] = cpu_hours_data[key]
				else:
					merged_results[key].update(cpu_hours_data[key])
	status.update_status_db(merged_results, server_type='uppmax')




