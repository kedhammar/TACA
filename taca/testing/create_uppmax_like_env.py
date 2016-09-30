""" Load and parse configuration file
"""

import glob
import logging
from taca.utils.config import CONFIG
import couchdb
import os
import datetime
from dateutil.relativedelta import relativedelta
import random

logger = logging.getLogger(__name__)


def setupServer(conf):
    url="http://{0}:{1}@{2}:{3}".format(conf['username'], conf['password'], conf['url'], conf['port'])
    return couchdb.Server(url)


def select_random_projects(projects_in, num_proj, application, projects_out):
    chosen_projects = 0
    projects_in_copy = projects_in
    application_not_in_other = ["WG re-seq"]
    while chosen_projects != num_proj and len(projects_in_copy) > 0:
        selected_proj = random.choice(projects_in_copy.keys())
        proj_value = projects_in_copy[selected_proj]
        del projects_in_copy[selected_proj] # remove so I am sure list decreases and the condiction of while is respected
        if application == "other":
            #in this case everything expcept
            if proj_value["application"] not in application_not_in_other:
                #I select this one
                projects_out.append(selected_proj)
        elif application == proj_value["application"]:
            #I select this one
            projects_out.append(selected_proj)




def create(projects, ngi_config):
    #connect to statusdb
    couch_info = CONFIG.get('statusdb')
    if couch_info is None:
        logger.error("No statusdb field in taca configuration file")
        return 1
    if "dev" not in couch_info["url"]:
        logger.error("url for status db is {}, but dev must be specified in this case".format(couch_info["url"]))
    couch=setupServer(couch_info)
    # connect to db and to view
    projectsDB = couch["projects"]
    project_summary = projectsDB.view("project/summary")
    projects_closed_more_than_three_months = {}
    projects_closed_more_than_one_month_less_than_three = {}
    projects_closed_less_than_one_month    = {}
    projects_opened = {}
    current_date =  datetime.datetime.today()
    date_limit_one_year = current_date - relativedelta(months=6) #yes yes I know.. but in this way i am sure all data in in xflocell_db
    date_limit_one_month = current_date - relativedelta(months=1)
    date_limit_three_month = current_date - relativedelta(months=3)
    for row in project_summary:
        project_id = row["key"][1]
        project_status = row["key"][0]
        if "application" not in row["value"]:
            continue
        application = row["value"]["application"]
        if project_status == "closed":
            if "close_date" in row["value"]:
                close_date = datetime.datetime.strptime(row["value"]["close_date"], '%Y-%m-%d')
                if close_date > date_limit_one_year: #if the project has been closed after the date limit
                    if close_date >= date_limit_one_month:
                        projects_closed_less_than_one_month[project_id] = {"project_name": row["value"]["project_name"],
                                                                            "application": application}
                    elif close_date < date_limit_one_month and close_date >= date_limit_three_month:
                        projects_closed_more_than_one_month_less_than_three[project_id] = {"project_name": row["value"]["project_name"],
                                                                            "application": application}
                    elif close_date < date_limit_three_month:
                        projects_closed_more_than_three_months[project_id] = {"project_name": row["value"]["project_name"],
                                                                            "application": application}
        elif project_status == "open":
            if "lanes_sequenced" in row["value"] and row["value"]["lanes_sequenced"] > 0:
                projects_opened[project_id] =  {"project_name": row["value"]["project_name"],
                                            "application": application}
        else:
            print "status {}".format(project_status)
    ##now I can parse the x_flowcell db to check what I can and cannot use
    ##it is less than one year we are using the flowcell_db so old projects might be not present
    whole_genome_projects = int(2*projects/3)
    import pdb
    pdb.set_trace()
    projects_to_reproduce = []
    select_random_projects(projects_closed_more_than_three_months, whole_genome_projects/4, "WG re-seq", projects_to_reproduce)
    select_random_projects(projects_closed_more_than_one_month_less_than_three, whole_genome_projects/4, "WG re-seq", projects_to_reproduce)
    select_random_projects(projects_closed_less_than_one_month,whole_genome_projects/4, "WG re-seq", projects_to_reproduce)
    select_random_projects(projects_opened, whole_genome_projects/4, "WG re-seq", projects_to_reproduce)

    other_projects = int(projects/3)
    select_random_projects(projects_closed_more_than_three_months, other_projects/4, "other", projects_to_reproduce)
    select_random_projects(projects_closed_more_than_one_month_less_than_three, other_projects/4, "other", projects_to_reproduce)
    select_random_projects(projects_closed_less_than_one_month, other_projects/4, "other", projects_to_reproduce)
    select_random_projects(projects_opened, other_projects/4, "other", projects_to_reproduce)

    ### At this point I scan over x_flowcell and reproduce FCs
    


    # fetch all FCs that contain this project --> discuss with denis, I can emit as many times as I want per doc, and later filter them in python
    #http://tools-dev.scilifelab.se:5984/_utils/database.html?x_flowcells/_design/unsafe/_view/projects_lanes
    # recreate the FCs (netirely, each fastq.gz is a real fastq files, the same)
    # print commands to be run (ngi_pipeline) to organise the projects
    # create the config_file for ngi_pipeline, most likely needed to export some varianble...
    # if ngi_pipeline available, organise stuff
    # if previous step worked, then create the ANALYSIS folder



