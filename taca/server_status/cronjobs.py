import logging
import subprocess
from crontab import CronTab
import platform
import getpass
import datetime

import couchdb

from taca.utils.config import CONFIG

def update_cronjob_db():
    logging.info('Getting list of cronjobs')
    crontab = CronTab()
    server = platform.node().split('.')[0]
    user = getpass.getuser()
    timestamp = datetime.datetime.now()
    result = []
    # parse results
    for job in crontab.crons:
        result.append({'Command': job.command,
                       'Comment': job.comment,
                       'Enabled': job.enabled,
                       'Minute': str(job.minutes),
                       'Hour': str(job.hours),
                       'Day of month' : str(job.month),
                       'Month': str(job.month),
                       'Day of week': str(job.day)})
    print result
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
        for row in view[server]:
            doc = crontab_db.get(row.value)
            doc['users'][user] = result
            doc['Last updated'] = str(timestamp)
            logging.info('Updating the document')
            try:
                crontab_db.save(doc)
            except Exception, e:
                logging.error(e.message)
            else:
                logging.info('Document {} has been successfully updated'.format(row.value))
