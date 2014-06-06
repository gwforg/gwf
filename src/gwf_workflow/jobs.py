"""
Module for keeping track of job status for jobs in the grid queue.
"""

import shelve
import subprocess
import os.path


def get_job_status(job_id):
    try:
        stat = subprocess.Popen(['qstat', '-f', job_id], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in stat.stdout:
            line = line.strip()
            if line.startswith('job_state'):
                status = line.split()[2]
                return status
    except:
        return False

def make_db_file_name(workflow_directory):
    return os.path.join(workflow_directory, '.jobs')

class JobsDatabase(object):
    """Wraps a jobs database for a given working directory."""

    def __init__(self, db_file):
        self.db = shelve.open(db_file)
        self.status_db = {}
        self._read_and_update_status()

    def _read_and_update_status(self):
        for job_name in self.db:
            job_id = self.db[job_name]
            job_status = get_job_status(job_id)
            print job_name
            if job_status is None:
                # It is no longer in the queue so it shouldn't be in the jobs db
                print 'deleting', job_name
                del self.db[job_name]
                self.db.sync()
            else:
                print job_name, 'has status', job_status
                self.status_db[job_name] = job_status
        self.db.close()

    def set_job_id(self, target_name, job_id):
        self.db[target_name] = job_id
        self.status_db[target_name] = 'Q' # It starts out as a queued object...

    def in_queue(self, job_name):
        #return job_name in self.status_db
        return self.db.has_key(job_name)

    def get_job_status(self, job_name):
        if job_name in self.status_db:
            return self.status_db[job_name]
        else:
            return None


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

# Global data base access for a running workflow...
JOBS_QUEUE = JobsDBCollection()
