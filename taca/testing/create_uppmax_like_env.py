""" Load and parse configuration file
"""

import glob
import logging
import os

logger = logging.getLogger(__name__)

def create(projects):
    #connect to statusdb
    #download from there the number of reqeusted projects
    # only projects open less than one year ago --> need a view in projects db
    #  - 2/3 Whole Human Genome
    #    - 1/4 closed more than 3 months ago
    #    - 1/4 closed more than 3 months ago and less than 1 months
    #    - 1/4 closed less than 1 month ago
    #    - 1/4 not closed
    #  - 1/3 others
    #    - 1/3 closed more than 3 months ago
    #    - 1/3 closed less than 1 months ago
    #    - 1/3 not closed
    # fetch all FCs that contain this project --> discuss with denis, I can emit as many times as I want per doc, and later filter them in python
    #http://tools-dev.scilifelab.se:5984/_utils/database.html?x_flowcells/_design/unsafe/_view/projects_lanes
    # recreate the FCs (netirely, each fastq.gz is a real fastq files, the same)
    # print commands to be run (ngi_pipeline) to organise the projects
    # create the config_file for ngi_pipeline, most likely needed to export some varianble...
    # if ngi_pipeline available, organise stuff
    # if previous step worked, then create the ANALYSIS folder
    print "hello world"
