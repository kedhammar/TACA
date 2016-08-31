import logging
import subprocess
from crontab import CronTab
import platform
import getpass
import datetime

import couchdb

from taca.utils.config import CONFIG

def _parse_crontab():
    result = {}
    list_of_users = CONFIG.get('server_status', {}).get('crontab_users', [])
    if not list_of_users:
        logging.warning('CONFIG.server_status does not have any crontab_users')
        current_user = getpass.getuser()
        list_of_users = [current_user]
        logging.warning('Using default user: {}'.format(current_user))
    for user in list_of_users:
        logging.info('Getting list of cronjobs for user: {}'.format(user))
        try:
            crontab = CronTab(user=user)
        except Exception, e:
            logging.error('Cannot get a crontab for user: {}'.format(user))
            logging.error(e.message)
        else:
            result[user] = []
            for job in crontab.crons:
                result[user].append({'Command': job.command,
                               'Comment': job.comment,
                               'Enabled': job.enabled,
                               'Minute': str(job.minutes),
                               'Hour': str(job.hours),
                               'Day of month' : str(job.month),
                               'Month': str(job.month),
                               'Day of week': str(job.day)})
    return result


def update_cronjob_db():
    server = platform.node().split('.')[0]
    timestamp = datetime.datetime.now()
    # parse results
    result = _parse_crontab()
    # connect to db
    url = "http://{username}:{password}@{url}:{port}".format(
            url=CONFIG.get('statusdb', {}).get('url'),
            username=CONFIG.get('statusdb', {}).get('username'),
            password=CONFIG.get('statusdb', {}).get('password'),
            port=CONFIG.get('statusdb', {}).get('port'))
    logging.info('Connecting to database: {}'.format(CONFIG.get('statusdb', {}).get('url')))
    try:
        couch = couchdb.Server(url)
    except Exception, e:
        logging.error(e.message)
    else:
        # update document
        crontab_db = couch['cronjobs']
        view = crontab_db.view('server/alias')
        # to be safe
        doc = {}
        # create doc if not exist
        if not view[server].rows:
            logging.info('Creating a document')
            doc = {
                'users': {user: result},
                'Last updated': str(timestamp),
                'server': server,
            }
        # else: get existing doc
        for row in view[server]:
            logging.info('Updating the document')
            doc = crontab_db.get(row.value)
            doc['users'] = result
            doc['Last updated'] = str(timestamp)
        if doc:
            try:
                crontab_db.save(doc)
            except Exception, e:
                logging.error(e.message)
            else:
                logging.info('{} has been successfully updated'.format(server))
        else:
            logging.warning('Document has not been created/updated')

