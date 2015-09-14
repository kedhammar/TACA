

import subprocess
import sys

import json
import gspread
from oauth2client.client import SignedJwtAssertionCredentials as GCredentials
import os
import logging

from config import SERVERS, USER, G_CREDENTIALS, G_SCOPE, G_SHEET, G_SHEET_MAP

def get_disk_space():
	result = {}
	for server_url in SERVERS.keys():
		# get server name, e.g.: milou.uppmax.uu.se -> milou
		server_name = server_url.split('.')[0]
		# insert path in command
		command = SERVERS[server_url]

		logging.debug(server_url)
		logging.debug(command)
		print server_url
		print command

		if server_url == "localhost":
			ssh = None
			for subcommand in command.split('|'):
				logging.debug(subcommand)
				stdin = ssh.stdout if ssh is not None else subprocess.PIPE
				ssh = subprocess.Popen(subcommand, stdin=stdin)
		else:
			# connect via ssh to server and execute the command
			ssh = subprocess.Popen(['ssh', '-t', '%s@%s' %(USER, server_url), command],
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE)
		# recieve how many percent of space left
		output = ssh.stdout.read()
		logging.debug(output)
		print output

		# from str to int
		used_space = int(output.strip().replace('%', ''))
		available_space = 100 - used_space
		result[server_name] = str(available_space) + '%'
	return result


def update_google_docs(data):
	# open json file
	json_key = json.load(open(G_CREDENTIALS))
	# get credentials from the file and authorize
	credentials = GCredentials(json_key['client_email'], json_key['private_key'], G_SCOPE)
	gc = gspread.authorize(credentials)
	# open the file
	# IMPORTANT: file must be shared with email listed in credentials
	sheet = gc.open(G_SHEET)

	# choose worksheet from the doc
	worksheet = sheet.get_worksheet(1)

	# update cell
	for key in data:
		cell = G_SHEET_MAP.get(key) # key = server name
		value = data.get(key)		# value = available space
		worksheet.update_acell(cell, value)