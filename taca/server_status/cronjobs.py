import logging
import platform
import getpass
import datetime

from crontab import CronTab
from taca.utils import statusdb
from taca.utils.config import CONFIG

def _parse_crontab():
    result = {}
    user = getpass.getuser()
    logging.info('Getting crontab for user {}'.format(user))
    try:
        crontab = CronTab(user=user)
    except Exception, e:
        logging.error('Cannot get a crontab for user: {}'.format(user))
        logging.error(e.message)
    else:
        result[user] = []
        for job in crontab.crons:
            # this is for special syntax like @monthly or @reboot
            special_syntax = str(job).split()[0] if str(job).startswith('@') else ''
            result[user].append({'Command': job.command,
                           'Comment': job.comment,
                           'Enabled': job.enabled,
                           'Minute': str(job.minutes),
                           'Hour': str(job.hours),
                           'Day of month' : str(job.dom),
                           'Month': str(job.month),
                           'Day of week': str(job.dow),
                           'Special syntax': special_syntax})
    return result


def update_cronjob_db():
    server = platform.node().split('.')[0]
    timestamp = datetime.datetime.now()
    # parse results
    result = _parse_crontab()
    # connect to db
    statusdb_conf = CONFIG.get('statusdb')
    logging.info('Connecting to database: {}'.format(CONFIG.get('statusdb', {}).get('url')))
    try:
        couch_connection = statusdb.StatusdbSession(statusdb_conf).connection
    except Exception, e:
        logging.error(e.message)
    else:
        # update document
        crontab_db = couch_connection['cronjobs']
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

