from __future__ import absolute_import
from __future__ import print_function

import os.path
import subprocess

import gwf

from gwf.backends.base import Backend


class LocalBackend(Backend):
    """Backend functionality for local execution."""

    def __init__(self):
        self.next_job_id = 0

    def get_state_of_jobs(self, job_ids):
        fake_table = dict(zip(job_ids, ["?"] * len(job_ids)))
        return fake_table

    def write_script_header(self, f, options):
        pass

    def write_script_variables(self, f):
        print('export GWF_JOBID={}'.format(self.next_job_id), file=f)

    def _build_submit_command(self, target, script_name, dependent_ids):
        command = ["bash", script_name]
        return command

    def submit_command(self, target, script_name, dependents_ids):
        log_dir = os.path.join(gwf.WORKING_DIR, 'gwf_log')
        command = self._build_submit_command(target, script_name, dependents_ids)
        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(os.path.join(log_dir, target.name + '.stdout'), "w") as outfile:
            print(qsub.stdout.read(), file=outfile)
        with open(os.path.join(log_dir, target.name + '.stderr'), "w") as outfile:
            print(qsub.stderr.read(), file=outfile)
        job_id = self.next_job_id
        self.next_job_id += 1
        return job_id

    def build_cancel_command(self, job_ids):
        pass
