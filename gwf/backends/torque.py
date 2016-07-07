from __future__ import absolute_import
from __future__ import print_function

import subprocess
import os
import os.path

import gwf

from gwf.backends.base import Backend
from gwf.helpers import make_absolute_path, escape_file_name


def _mkdir_if_not_exist(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


class TorqueBackend(Backend):
    """Backend functionality for torque."""

    def __init__(self):
        self.script_dir = make_absolute_path(gwf.WORKING_DIR, '.scripts')

    def get_state_of_jobs(self, job_ids):
        result = dict((job_id, None) for job_id in job_ids)
        for job_id in job_ids:
            try:
                stat = subprocess.Popen(['qstat', '-f', job_id],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
                for line in stat.stdout:
                    line = line.strip()
                    if line.startswith('job_state'):
                        status = line.split()[2]
                        result[job_id] = status
            except:
                pass
        return result

    def write_script_header(self, f, options):
        print('#PBS -l nodes={}:ppn={}'.format(
            options['nodes'], options['cores']), file=f)
        print('#PBS -l mem={}'.format(options['memory']), file=f)
        print('#PBS -l walltime={}'.format(options['walltime']), file=f)
        if 'queue' in options:
            print('#PBS -q {}'.format(options['queue']), file=f)
        if 'account' in options:
            print('#PBS -A {}'.format(options['account']), file=f)

    def write_script_variables(self, f):
        print('export GWF_JOBID=$PBS_JOBID', file=f)

    def _build_submit_command(self, target, script_name, dependents_ids):
        self.script_name = make_absolute_path(self.script_dir, escape_file_name(target.name))

        log_dir = os.path.join(gwf.WORKING_DIR, 'gwf_log')
        _mkdir_if_not_exist(log_dir)

        command = ['qsub', '-N', target.name,
                   '-o', os.path.join(log_dir, target.name + '.stdout'),
                   '-e', os.path.join(log_dir, target.name + '.stderr'),
                   ]
        if len(dependents_ids) > 0:
            command.append('-W')

            ids = ':'.join(dependents_ids)
            command.append('depend=afterok:{}'.format(ids))
        command.append(script_name)
        return command

    def submit_command(self, target, script_name, dependents_ids):
        command = self._build_submit_command(target, script_name, dependents_ids)
        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        job_id = qsub.stdout.read().strip()
        return job_id

    def build_cancel_command(self, job_ids):
        return ['qdel'] + map(str, job_ids)
