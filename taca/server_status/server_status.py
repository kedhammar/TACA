import subprocess
import json
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
        data can be from nases or from uppmax
        server_type should be either 'uppmax' or 'nas'
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
        server['name'] = key # key is nas url or uppmax project
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

def get_uppmax_quotas():
    current_time = datetime.datetime.now()
    try:
        uq = subprocess.Popen(["/sw/uppmax/bin/uquota", "-q"], stdout=subprocess.PIPE)
    except Exception, e:
        logging.error(e.message)
        raise e

    output = uq.communicate()[0]
    logging.info("Disk Usage:")
    logging.info(output)

    projects = output.split("\n/proj/")[1:]

    result = {}
    for proj in projects:
        project_dict = {"time": current_time.isoformat()}
        project = proj.strip("\n").split()
        project_dict["project"] = project[0]
        project_dict["usage (GB)"] = project[1]
        project_dict["quota limit (GB)"] = project[2]
        try:
            project_dict["quota_decrease"] = project[3]
        except:
            pass

        result[project[0]] = project_dict
    return result

def get_uppmax_cpu_hours():
    current_time = datetime.datetime.now()
    try:
        # script that runs on uppmax
        uq = subprocess.Popen(["/sw/uppmax/bin/projinfo", '-q'], stdout=subprocess.PIPE)
    except Exception, e:
        logging.error(e.message)
        raise e

    # output is lines with the format: project_id  cpu_usage  cpu_limit
    output = uq.communicate()[0]

    logging.info("CPU Hours Usage:")
    logging.info(output)
    result = {}
    # parsing output
    for proj in output.strip().split('\n'):
        project_dict = {"time": current_time.isoformat()}

        try: # split line into a list
            project_id, cpu_hours, cpu_limit = proj.split()
        # sometimes it returns empty strings or something strange
        except Exception, e:
            logging.error(e.message)
        else:
            project_dict["project"] = project_id
            project_dict["cpu hours"] = cpu_hours
            project_dict["cpu limit"] = cpu_limit

            quotas = _get_uppmax_cpu_quotas(project_id)
            if quotas:
                project_dict["cpu_quota_decrease"] = quotas

            result[project_id] = project_dict

    return result

def _get_uppmax_cpu_quotas(project_id):
    today = datetime.date.today()
    result = []
    try:
        logging.info("CPU Quotas Decrease, project {}:".format(project_id))
        cpu_quotas = subprocess.Popen("grep -A 10 {} /sw/uppmax/etc/projects".format(project_id).split())
    except Exception, e:
        logging.error(e.message)
    else:
        output = cpu_quotas.communicate()[0]
        for line in output:
            # line should be: "Shortgrants:  milou=7000	20150421"
            if 'Shortgrants' in line:
                logging.info(line)
                try:
                    shortgrants, cpu_hours, date_string = line.split()
                except Exception, e:
                    logging.error(e.message)
                else:
                    quota_date = datetime.datetime.strptime(date_string, "%Y%m%d").date()
                    if quota_date > today:
                        # cpu_hours will be 'milou=20000' or 'irma=20000'
                        cpu_hours = cpu_hours.split('='[-1])
                        result.append({'quota': cpu_hours, 'date': quota_date})
    return result
