"""
Script to move job status information from the old .jobs directory to a database.

It doesn't do this recursively based on a workflow but only a single .jobs/ directory, so be careful
when porting a running workflow!
"""

import os.path
from os import listdir
import sys
import shelve

if len(sys.argv) > 1:
    workflow_dir = sys.argv[1]
else:
    workflow_dir = '.'

jobs_dir = os.path.join(workflow_dir, '.jobs')
jobs = listdir(jobs_dir)
jobs_db = shelve.open(os.path.join(workflow_dir, '.jobs'))

for job_name in jobs:
    job_file = os.path.join(jobs_dir, job_name)
    print 'job file for', job_name, job_file

    try:
        print open(job_file).read().strip()
        job_id = int(open(job_file).read().strip())
        print 'job ID', job_id
        jobs_db[job_name] = job_id
    except:
        pass

jobs_db.close()

jobs_db = shelve.open(os.path.join(workflow_dir, '.jobs'))
for job_name in jobs_db:
    print job_name, jobs_db[job_name]