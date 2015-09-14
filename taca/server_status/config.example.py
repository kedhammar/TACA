import os

# key: server url, value: command
# df -h {path} - shows information about home directory as a table
# awk - takes the prelast column of the table (how many percent)
# sed - takes the second line of the output
SERVERS = {
		'example.com': "df -h /path/to/dir | awk '{print $(NF-1)}' | sed -n '$p'" ,
		'example2.com': "df -h /another/path/to/another/dir | awk '{print $(NF-1)}' | sed -n '$p'",
}

# for access to server without the password
# public key must be generated on client
# and copied to ~/.ssh/authorized_keys of the server
USER = "username"

# path to credentials file, must be configured on console.developers.google.com
G_CREDENTIALS = os.path.join(os.path.dirname(__file__), 'gdocs_credentials.json')
# name of the google doc: file must be shared with email listed in credentials file
G_SHEET = 'the google sheet'
# whatever default value - do not change
G_SCOPE = ['https://spreadsheets.google.com/feeds']
# addresses of the cells for each server -> where to insert values
G_SHEET_MAP = {	'server1': 'G22',
				'server2': 'C22'}