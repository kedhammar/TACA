import subprocess
import json
import gspread
from oauth2client.client import SignedJwtAssertionCredentials as GCredentials
from taca.utils.config import CONFIG

def get_disk_space():
	result = {}
	config = CONFIG['server_status']
	servers = config.get('servers', [])
	for server_url in servers.keys():
		# get path of disk
		path = servers[server_url]

		# get command
		command = "{command} {path}".format(command=config['command'], path=path)

		# if localhost, don't connect to ssh
		if server_url == "localhost":
			proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		else:
			# connect via ssh to server and execute the command
			proc = subprocess.Popen(['ssh', '-t', '%s@%s' %(config['user'], server_url), command],
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE)
		# read output
		output = proc.stdout.read()
		# parse output
		output = __parse_output(output)

		try:
			# remove % symbol and convert string to int
			used_space = int(output.strip().replace('%', ''))
			# how many percent available
			available_space = 100 - used_space
			result[server_url] = "{value}%".format(value=available_space)
		except:
			# sometimes it fails for whatever reason as Popen returns not what it is supposed to
			result[server_url] = 'NaN'
	return result

def __parse_output(output):
	# command = df -h /home
	# output = Filesystem            Size  Used Avail Use% Mounted on
	# /dev/mapper/VGStor-lv_illumina
    #                   24T   12T   13T  49% /srv/illumina

	output = output.strip() # remove \n in the end
	output = output.split('\n')[-1] # split by lines and take the last line
	output = output.strip() # remove spaces
	output = output.split() # split line by space symbols

	# output = ['24T', '12T', '13T', '49%', '/srv/illumina']
	for item in output: # select the item containing '%' symbol
		if '%' in item:
			return item

	return 'NaN' # if no '%' in output, return NaN

def update_google_docs(data, credentials_file):
	config = CONFIG['server_status']
	# open json file
	json_key = json.load(open(credentials_file))

	# get credentials from the file and authorize
	credentials = GCredentials(json_key['client_email'], json_key['private_key'], config['g_scope'])
	gc = gspread.authorize(credentials)
	# open google sheet
	# IMPORTANT: file must be shared with email listed in credentials
	sheet = gc.open(config['g_sheet'])

	# choose worksheet from the doc
	worksheet = sheet.get_worksheet(1)

	# update cell
	for key in data:
		cell = config['g_sheet_map'].get(key) # key = server name
		value = data.get(key)		# value = available space
		worksheet.update_acell(cell, value)