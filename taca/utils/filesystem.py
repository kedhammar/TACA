""" Filesystem utilities
"""
import contextlib
import os
import re
import shutil
from subprocess import check_call, CalledProcessError, Popen, PIPE

RUN_RE = '^\d{6}_[a-zA-Z\d\-]+_\d{4}_[AB0][A-Z\d\-]+$'
PROJECT_RE = '[a-zA-Z]+\.[a-zA-Z]+_\d{2}_\d{2}'

@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.
    """
    cur_dir = os.getcwd()
    # This is weird behavior. I'm removing and and we'll see if anything breaks.
    #safe_makedir(new_dir)
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)

def create_folder(target_folder):
    """ Ensure that a folder exists and create it if it doesn't, including any
        parent folders, as necessary.

        :param target_folder: the target folder
        :returns: True if the folder exists or was created, False if the folder
        does not exists and could not be created
    """
    try:
        os.makedirs(target_folder)
    except OSError as e:
        pass
    return os.path.exists(target_folder)

def touch(file):
    open(file, "w").close()

def do_symlink(src_file, dst_file):
    link_f = os.symlink
    if not os.path.isfile(dst_file):
        link_f(os.path.realpath(src_file), dst_file)

def do_copy(src_path, dst_path):
    # copies folder structure and files (recursively)
    # if symlinks, will copy content, not the links
    # dst_path will be created, it must NOT exist
    shutil.copytree(src_path, dst_path)
