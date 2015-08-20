import subprocess
import sys

import json
import gspread
from oauth2client.client import SignedJwtAssertionCredentials as GCredentials
import os


# for access to server without entering the password
# public key must be generated on client 
# and copied to ~/.ssh/authorized_keys of the server 
SERVERS = ["milou.uppmax.uu.se", "nestor.uppmax.uu.se"]
USER = "ekat"
# df -h /home - shows information about home directory as a table
# awk - takes the prelast column of the table (how many percent) 
# sed - takes the second line of the output 
COMMAND = """df -h /home | awk '{print $(NF-1)}' | sed -n 2p"""

# credentials file, must be configure on console.developers.google.com
G_CREDENTIALS = os.path.join(os.getcwd(), 'gdocs_credentials.json')
# name of the google doc
G_SHEET = 'COPY'
# whatever default value
G_SCOPE = ['https://spreadsheets.google.com/feeds']
# to store address of the cell for each server
G_SHEET_MAP = {	'milou': 'G22',
				'nestor': 'E22'}

def get_disk_space():
	result = {}
	for server in SERVERS:
		# connect via ssh to server and execute the command
		ssh = subprocess.Popen(['ssh', '-t', '%s@%s' %(USER, server), COMMAND],
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE)
		# recieve how many percent of space left
		output = ssh.stdout.read()
		# get server name, e.g.: milou.uppmax.uu.se -> milou
		server_name = server.split('.')[0] 
		# from str to int
		used_space = int(output.replace('%', ''))
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
		cell = G_SHEET_MAP[key] # key = server name
		value = data[key]		# value = available space
		worksheet.update_acell(cell, value)
	print 'Document has been succesfully updated'
