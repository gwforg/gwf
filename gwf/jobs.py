"""
Module for keeping track of job status for jobs in the grid queue.
"""

import os.path
import shelve
from itertools import izip

import gwf


def make_db_file_name(workflow_directory):
    return os.path.join(workflow_directory, '.jobs')


class JobsDatabase(object):
    """Wraps a jobs database for a given working directory."""

    def __init__(self, db_file):
        self.db = shelve.open(db_file)
        self.status_db = {}
        self._read_and_update_status()

    def _read_and_update_status(self):
        job_names = list(self.db.keys())
        job_ids = list(self.db[name] for name in job_names)
        current_state = gwf.backends.BACKEND.get_state_of_jobs(job_ids)
        for job_name, job_id in izip(job_names, job_ids):
            job_id = self.db[job_name]
            job_status = current_state[job_id]
            if job_status not in ('Q', 'R', 'H'):
                # It is no longer in the queue so it shouldn't be in the jobs db
                del self.db[job_name]
                self.db.sync()
            else:
                self.status_db[job_name] = job_status

    def set_job_id(self, target_name, job_id):
        self.db[target_name] = job_id
        self.status_db[target_name] = 'Q'  # It starts out as a queued object...

    def in_queue(self, job_name):
        return job_name in self.status_db

    def get_job_id(self, job_name):
        return self.db.get(job_name)

    def get_job_status(self, job_name):
        return self.status_db.get(job_name)

    def close(self):
        self.db.close()


class JobsDBCollection(object):
    """Collects databases for several working directories."""

    def __init__(self):
        self.databases = {}

    def get_database(self, workflow_directory):
        """
        :param workflow_directory: Directory of the workflow the job was specified in.
        :return: The jobs database associated to the workflow directory
        :rtype: gwf_workflow.jobs.JobsDatabase
        """
        if workflow_directory not in self.databases:
            self.databases[workflow_directory] = JobsDatabase(make_db_file_name(workflow_directory))
        return self.databases[workflow_directory]

    def close(self):
        for db in self.databases.values():
            db.close()


# Global data base access for a running workflow...
JOBS_QUEUE = JobsDBCollection()

# Necessary to close all databases at exit...
import atexit


@atexit.register
def close_databases():
    JOBS_QUEUE.close()
