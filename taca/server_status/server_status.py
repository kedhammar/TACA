import subprocess
import logging
import couchdb
import datetime

from taca.utils.config import CONFIG

def get_nases_disk_space():
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
            proc = subprocess.Popen(['ssh', '-t', '{}@{}'.format(config['user'], server_url), command],
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
        output = proc.stdout.read()
        output = _parse_output(output)
        result[server_url] = output
    return result

def _parse_output(output): # for nases
    # command = df -h /home
    # output = Filesystem            Size  Used Avail Use% Mounted on
    # /dev/mapper/VGStor-lv_illumina
    #                   24T   12T   13T  49% /srv/illumina

    output = output.strip()
    output = output.split()
    try:
        mounted_on = output[-1]
        used_percentage = output[-2]
        space_available = output[-3]
        space_used = output[-4]
        disk_size = output[-5]
        filesystem = output[-6]

        available_percentage = str(100 - int(used_percentage.replace('%',''))) + '%'

        result = {
            'disk_size': disk_size,
            'space_used': space_used,
            'space_available': space_available,
            'used_percentage': used_percentage,
            'available_percentage': available_percentage,
            'mounted_on': mounted_on,
            'filesystem': filesystem
        }
    except:
        # sometimes it fails for whatever reason as Popen returns not what it is supposed to
        result = {
            'disk_size': 'NaN',
            'space_used': 'NaN',
            'space_available': 'NaN',
            'used_percentage': 'NaN',
            'available_percentage': 'NaN',
            'mounted_on': 'NaN',
            'filesystem': 'NaN'
        }
        logging.error("Can't parse the output: {}".format(output))

    return result

def update_status_db(data, server_type=None):
    """ Pushed the data to status db,
        data can be from nases
        server_type should be 'nas'
    """
    db_config = CONFIG.get('statusdb')
    if db_config is None:
        logging.error("'statusdb' must be present in the config file!")
        raise RuntimeError("'statusdb' must be present in the config file!")

    server = "http://{username}:{password}@{url}:{port}".format(
        url=db_config['url'],
        username=db_config['username'],
        password=db_config['password'],
        port=db_config['port'])
    try:
        couch = couchdb.Server(server)
    except Exception, e:
        logging.error(e.message)
        raise

    db = couch['server_status']
    logging.info('Connection established')
    for key in data.keys(): # data is dict of dicts
        server = data[key] # data[key] is dictionary (the command output)
        server['name'] = key # key is nas url
        # datetime.datetime(2015, 11, 18, 9, 54, 33, 473189) is not JSON serializable
        server['time'] = datetime.datetime.now().isoformat()
        server['server_type'] = server_type or 'unknown'

        try:
            db.save(server)
        except Exception, e:
            logging.error(e.message)
            raise
        else:
            logging.info('{}: Server status has been updated'.format(key))
