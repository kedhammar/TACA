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
    user = getpass.getuser()
    logging.info('Getting crontab for user {}'.format(user))
    crontab = None
    try:
        crontab = CronTab(user=user)
    except Exception, e:
        logging.warning('Cannot get a crontab for user: {}'.format(user))
        logging.warning(e.message)
        try:
            logging.info('Getting a crontab for default user')
            crontab = CronTab()
        except Exception, e:
            logging.error('Cannot get a crontab for default user')
            logging.error(e.message)
    if crontab:
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
    print result
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
                'users': {user: cronjobs for user, cronjobs in result.items()},
                'Last updated': str(timestamp),
                'server': server,
            }
        # else: get existing doc
        for row in view[server]:
            logging.info('Updating the document')
            doc = crontab_db.get(row.value)
            doc['users'].update(result)
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

