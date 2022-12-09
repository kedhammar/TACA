"""Classes for handling connection to StatusDB."""

import couchdb
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

class StatusdbSession(object):
    """Wrapper class for couchdb."""
    def __init__(self, config, db=None):
        user = config.get('username')
        password = config.get('password')
        url = config.get('url')
        url_string = 'https://{}:{}@{}'.format(user, password, url)
        display_url_string = 'https://{}:{}@{}'.format(user, '*********', url)
        self.connection = couchdb.Server(url=url_string)
        if not self.connection:
            raise Exception('Couchdb connection failed for url {}'.format(display_url_string))
        if db:
            self.db_connection = self.connection[db]

    def get_entry(self, name, use_id_view=False):
        """Retrieve entry from a given db for a given name.

        :param name: unique name identifier (primary key, not the uuid)
        """
        if use_id_view:
            view = self.id_view
        else:
            view = self.name_view
        if not view.get(name, None):
            return None
        return self.db.get(view.get(name))

    def save_db_doc(self, doc, db=None):
        try:
            db = db or self.db
            db.save(doc)
        except Exception as e:
            raise Exception('Failed saving document due to {}'.format(e))

    def get_project_flowcell(self, project_id, open_date='2015-01-01', date_format='%Y-%m-%d'):
        """From information available in flowcell db connection,
        collect the flowcell this project was sequenced.

        :param project_id: NGI project ID to get the flowcells
        :param open_date: Open date of project to skip the check for all flowcells
        :param date_format: The format of specified open_date
        """
        try:
            open_date = datetime.strptime(open_date, date_format)
        except:
            open_date = datetime.strptime('2015-01-01', '%Y-%m-%d')

        project_flowcells = {}
        date_sorted_fcs = sorted(list(self.proj_list.keys()), key=lambda k: datetime.strptime(k.split('_')[0], '%y%m%d'), reverse=True)
        for fc in date_sorted_fcs:
            fc_date, fc_name = fc.split('_')
            if datetime.strptime(fc_date,'%y%m%d') < open_date:
                break
            if project_id in self.proj_list[fc] and fc_name not in project_flowcells.keys():
                project_flowcells[fc_name] = {'name':fc_name,'run_name':fc, 'date':fc_date, 'db':self.db.name}
        return project_flowcells

class ProjectSummaryConnection(StatusdbSession):
    def __init__(self, config, dbname='projects'):
        super(ProjectSummaryConnection, self).__init__(config)
        self.db = self.connection[dbname]
        self.name_view = {k.key: k.id for k in self.db.view('project/project_name', reduce=False)}
        self.id_view = {k.key: k.id for k in self.db.view('project/project_id', reduce=False)}


class FlowcellRunMetricsConnection(StatusdbSession):
    def __init__(self, config, dbname='flowcells'):
        super(FlowcellRunMetricsConnection, self).__init__(config)
        self.db = self.connection[dbname]
        self.name_view = {k.key:k.id for k in self.db.view('names/name', reduce=False)}
        self.proj_list = {k.key:k.value for k in self.db.view('names/project_ids_list', reduce=False) if k.key}


class X_FlowcellRunMetricsConnection(StatusdbSession):
    def __init__(self, config, dbname='x_flowcells'):
        super(X_FlowcellRunMetricsConnection, self).__init__(config)
        self.db = self.connection[dbname]
        self.name_view = {k.key:k.id for k in self.db.view('names/name', reduce=False)}
        self.proj_list = {k.key:k.value for k in self.db.view('names/project_ids_list', reduce=False) if k.key}

class NanoporeRunsConnection(StatusdbSession):
    
    def __init__(self, config, dbname='nanopore_runs'):
        super(NanoporeRunsConnection, self).__init__(config)
        self.db = self.connection[dbname]

    def update_db(self, run_path_str, dict2add = None):
        """ Use the run path string "experiment_dir/sample_dir/run_dir" string to find any 
        matching entries in CouchDB. 
        
        If such an entry is found, marked as "ongoing" and there is info to add
        --> update it with dict2add, containing the finished run info.
        
        If no such entry is found
        --> create it.
        """

        try:
            path2stat = self.db.view('info/run_status')
            matching_rows = path2stat[run_path_str].rows

            # If matching entry exists in db
            if len(matching_rows) == 1:
                
                # If there is a dict to add to the entry
                if dict2add:
                    # If the entry is an ongoing run
                    if doc["run_status"] == "ongoing":
                        # Fetch run document from database
                        doc_id = matching_rows[0].id
                        doc = self.db[doc_id]

                        # Add finished run information to document and change status
                        doc.update(dict2add)
                        doc["run_status"] = "finished"

                        # Overwrite the database entry
                        self.db[doc_id] = doc
                        logger.info(f"Run report .json appended to database entry {run_path_str}, id {doc_id}")
                    else:
                        logger.info(f"Database entry found for {run_path_str}, but run is not ongoing")
                else:
                    logger.info(f"Database entry found for {run_path_str}, but no new information to add")

            elif len(matching_rows) == 0:
                # Create ongoing run       
                new_doc = {"run_path": run_path_str, "run_status": "ongoing"}
                new_doc_id, new_doc_rev = self.db.save(new_doc)
                logger.info(f"New database entry created: {run_path_str}, id {new_doc_id}, rev {new_doc_rev}")
            
            # Multiple matching rows
            else:
                logger.warn(f'More than one database entry with run_path {run_path_str} found')
                return False
            
            return True
        
        except:
            return False


def update_doc(db, obj, over_write_db_entry=False):
    view = db.view('info/name')
    if len(view[obj['name']].rows) == 1:
        remote_doc = view[obj['name']].rows[0].value
        doc_id = remote_doc.pop('_id')
        doc_rev = remote_doc.pop('_rev')
        if remote_doc != obj:
            if not over_write_db_entry:
                obj = merge_dicts(obj, remote_doc)
            obj['_id'] = doc_id
            obj['_rev'] = doc_rev
            db[doc_id] = obj
            logger.info('Updating {}'.format(obj['name']))
    elif len(view[obj['name']].rows) == 0:
        db.save(obj)
        logger.info('Saving {}'.format(obj['name']))
    else:
        logger.warn('More than one row with name {} found'.format(obj['name']))


def merge_dicts(d1, d2):
    """Merge dictionary d2 into dictionary d1.
    If the same key is found, the one in d1 will be used.
    """
    for key in d2:
        if key in d1:
            if isinstance(d1[key], dict) and isinstance(d2[key], dict):
                merge(d1[key], d2[key])
            elif d1[key] == d2[key]:
                pass  # same leaf value
            else:
                logger.debug('Values for key {key} in d1 and d2 differ, '
                          'using the value of d1'.format(key=key))
        else:
            d1[key] = d2[key]
    return d1
