""" Miscellaneous or general-use methods
"""
import couchdb
import hashlib
import os
import smtplib
import subprocess
import sys
import glob

from datetime import datetime
from email.mime.text import MIMEText


def send_mail(subject, content, receiver):
    """ Sends an email

    :param str subject: Subject for the email
    :param str content: Content of the email
    :param str receiver: Address to send the email
    """
    if not receiver:
        raise SystemExit("No receiver was given to send mail")
    msg = MIMEText(content)
    msg['Subject'] = "TACA - {}".format(subject)
    msg['From'] = 'TACA@scilifelab.se'
    msg['to'] = receiver

    s = smtplib.SMTP('localhost')
    s.sendmail('TACA', [receiver], msg.as_string())
    s.quit()

def call_external_command(cl, with_log_files=False, prefix=None, log_dir=""):
    """ Executes an external command

    :param string cl: Command line to be executed (command + options and parameters)
    :param bool with_log_files: Create log files for stdout and stderr
    :param string prefix: the prefics to add to log file
    :param string log_dir: where to write the log file (to avoid problems with rights)
    """
    if type(cl) == str:
        cl = cl.split(' ')
    logFile = os.path.basename(cl[0])
    stdout = sys.stdout
    stderr = sys.stderr
    if with_log_files:
        if prefix:
            logFile = '{}_{}'.format(prefix, logFile)
        #create log dir if it dint exist in CWD
        if log_dir and not os.path.exists(log_dir):
            os.mkdir(log_dir)
        logFile = os.path.join(log_dir, logFile)
        stdout = open(logFile + '.out', 'wa')
        stderr = open(logFile + '.err', 'wa')
        started = "Started command {} on {}".format(' '.join(cl), datetime.now())
        stdout.write(started + '\n')
        stdout.write(''.join(['=']*len(cl)) + '\n')

    try:
        subprocess.check_call(cl, stdout=stdout, stderr=stderr)
    except subprocess.CalledProcessError, e:
        e.message = "The command {} failed.".format(' '.join(cl))
        raise e
    finally:
        if with_log_files:
            stdout.close()
            stderr.close()

def call_external_command_detached(cl, with_log_files=False, prefix=None):
    """ Executes an external command

        :param string cl: Command line to be executed (command + options and parameters)
        :param bool with_log_files: Create log files for stdout and stderr
        """
    if type(cl) == str:
        cl = cl.split(' ')
    command = os.path.basename(cl[0])
    stdout = sys.stdout
    stderr = sys.stderr

    if with_log_files:
        if prefix:
            command = '{}_{}'.format(prefix, command)
        stdout = open(command + '.out', 'wa')
        stderr = open(command + '.err', 'wa')
        started = "Started command {} on {}".format(' '.join(cl), datetime.now())
        stdout.write(started + '\n')
        stdout.write(''.join(['=']*len(cl)) + '\n')

    try:
        p_handle = subprocess.Popen(cl, stdout=stdout, stderr=stderr)
    except subprocess.CalledProcessError, e:
        e.message = "The command {} failed.".format(' '.join(cl))
        raise e
    finally:
        if with_log_files:
            stdout.close()
            stderr.close()
    return p_handle

def to_seconds(days=None, hours=None):
    """ Convert given day/hours to seconds and return

        :param int days: days to be converted
        :param int hours: hours to be converted
    """
    #only either days or hours should be speciefied, but atleast one should be specified
    if days and hours:
        raise SystemExit('Both "days" and "hours" were given, use only either of them')
    elif not days and not hours:
        raise SystemExit('provide either "days" or "hours"')
    elif days and not hours:
        # 1 day == 60*60*24 seconds --> 86400
        return 86400 * days
    elif hours and not days:
        # 1 hour == 60*60 seconds --> 3600
        return 3600 * hours

def hashfile(afile, hasher='sha1', blocksize=65536):
    """ Calculate the hash digest of a file with the specified algorithm and
        return it.

        This solution was adapted from http://stackoverflow.com/a/3431835

        :param string afile: the file to calculate the digest for
        :param string hasher: the hashing algorithm to be used, default is sha1
        :param int blocksize: the blocksize to use, default is 65536 bytes
        :returns: the hexadecimal hash digest or None if input was not a file
    """
    if not os.path.isfile(afile):
        return None
    hashobj = hashlib.new(hasher)
    with open(afile,'rb') as fh:
        buf = fh.read(blocksize)
        while len(buf) > 0:
            hashobj.update(buf)
            buf = fh.read(blocksize)
    return hashobj.hexdigest()

def query_yes_no(question, default="yes", force=False):
    """Ask a yes/no question via raw_input() and return their answer.
    "question" is a string that is presented to the user. "default"
    is the presumed answer if the user just hits <Enter>. It must be
    "yes" (the default), "no" or None (meaning an answer is required
    of the user). The force option simply sets the answer to default.
    The "answer" return value is one of "yes" or "no".

    :param question: the displayed question
    :param default: the default answer
    :param force: set answer to default
    :returns: yes or no
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        if not force:
            choice = raw_input().lower()
        else:
            choice = "yes"
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                                 "(or 'y' or 'n').\n")

def return_unique(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def run_is_demuxed(run, couch_info=None):
    """Check in StatusDB 'x_flowcells' database if the given run has an entry which means it was
    demultiplexed (as TACA only creates a document upon successfull demultiplexing)

    :param str run: run name
    :param dict couch_info: a dict with 'statusDB' info
    """
    # check if statusdb info is given
    if not couch_info:
        raise SystemExit("To check for demultiplexing is enabled in config file but no 'statusDB' info was given")
    run_terms = run.split('_')
    run_date = run_terms[0]
    run_fc = run_terms[-1]
    run_name = "{}_{}".format(run_date, run_fc)
    # connect to statusdb using info fectched from config file
    try:
        server = "http://{username}:{password}@{url}:{port}".format(**couch_info)
        couch = couchdb.Server(server)
        fc_db = couch[couch_info['db']]
        for fc in fc_db.view("names/name", reduce=False, descending=True):
            if fc.key != run_name:
                continue
            fc_doc = fc_db.get(fc.id)
            if not fc_doc or not fc_doc.get('illumina', {}).get('Demultiplex_Stats', {}):
                return False
            return True
    except Exception, e:
        raise e
