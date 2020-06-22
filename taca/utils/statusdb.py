"""Classes for handling connection to StatusDB."""

import couchdb


class StatusdbSession(object):
    """Wrapper class for couchdb."""
    def __init__(self, config, db=None):
        user = config.get('username')
        password = config.get('password')
        port = config.get('port')
        url = config.get('url')
        url_string = 'http://{}:{}@{}:{}'.format(user, password, url, port)
        display_url_string = 'http://{}:{}@{}:{}'.format(user, '*********', url, port)
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
            if self.log:
                self.log.warn('no entry "{}" in {}'.format(name, self.db))
            return None
        return self.db.get(view.get(name))


class ProjectSummaryConnection(StatusdbSession):
    def __init__(self, config, dbname='projects'):
        super(ProjectSummaryConnection, self).__init__(config)
        self.db = self.connection[dbname]
        self.name_view = {k.key: k.id for k in self.db.view('project/project_name', reduce=False)}
        self.id_view = {k.key: k.id for k in self.db.view('project/project_id', reduce=False)}
