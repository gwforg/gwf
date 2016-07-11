from __future__ import absolute_import
from __future__ import print_function

import os
import os.path
import subprocess

import gwf

from gwf.backends.base import Backend
from gwf.helpers import make_absolute_path, escape_file_name


def _mkdir_if_not_exist(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


class SlurmBackend(Backend):
    """Backend functionality for slurm."""

    def _write_script(self, target):
        """Write the code to a script that can be executed."""

        script_dir = make_absolute_path(self.working_dir, '.scripts')
        script_name = make_absolute_path(self.script_dir, escape_file_name(target.name))

        if not os.path.exists(script_dir):
            os.makedirs(script_dir)

        with open(script_name, 'w') as f:
            print("#!/bin/bash", file=f)

            print('#SBATCH -N {}'.format(target.options['nodes']), file=f)
            print('#SBATCH -c {}'.format(target.options['cores']), file=f)
            print('#SBATCH --mem={}'.format(target.options['memory']), file=f)
            print('#SBATCH -t {}'.format(target.options['walltime']), file=f)
            if 'queue' in target.options:
                print('#SBATCH -p {}'.format(target.options['queue']), file=f)
            if 'account' in target.options:
                print('#SBATCH -A {}'.format(target.options['account']), file=f)
            if 'constraint' in target.options:
                print('#SBATCH -C {}'.format(target.options['constraint']), file=f)
            if 'mail_type' in target.options:
                print('#SBATCH --mail-type={}'.format(target.options['mail_type']), file=f)
            if 'mail_user' in target.options:
                print('#SBATCH --mail-user={}'.format(target.options['mail_user']), file=f)

            print(file=f)

            print('# GWF generated code ...', file=f)
            print('cd %s' % target.working_dir, file=f)
            print('export GWF_JOBID=$SLURM_JOBID', file=f)
            print(f, "set -e", file=f)
            print(file=f)

            print('# Script from workflow', file=f)
            print(target.spec, file=f)

    def get_state_of_jobs(self, job_ids):
        result = dict((job_id, False) for job_id in job_ids)
        map_state = {  # see squeue man page under JOB STATE CODES
            'BF': '?',  # BOOT_FAIL
            'CA': '?',  # CANCELLED
            'CD': '?',  # COMPLETED
            'CF': 'R',  # CONFIGURING
            'CG': 'R',  # COMPLETING
            'F': '?',  # FAILED
            'NF': '?',  # NODE_FAIL
            'PD': 'Q',  # PENDING
            'PR': '?',  # PREEMPTED
            'R': 'R',  # RUNNING
            'S': 'R',  # SUSPENDED
            'TO': '?',  # TIMEOUT
            'SE': 'Q',  # SPECIAL_EXIT
        }
        try:
            stat = subprocess.Popen(['squeue',
                                     '--noheader',
                                     '--format=%i;%t;%E'],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
            for line in stat.stdout:
                job_id, state, depends = line.strip().split(';')
                simple_state = map_state[state]
                if simple_state == 'Q' and depends != '':
                    result[job_id] = 'H'
                else:
                    result[job_id] = simple_state
        except:
            pass
        return result

    def _build_submit_command(self, target, script_name, dependents_ids):
        log_dir = os.path.join(target.working_dir, 'gwf_log')
        _mkdir_if_not_exist(log_dir)

        command = ['sbatch', '-J', target.name, '--parsable',
                   '-o', os.path.join(log_dir, target.name + '.%j.stdout'),
                   '-e', os.path.join(log_dir, target.name + '.%j.stderr'),
                   ]
        if len(dependents_ids) > 0:
            command.append('-d')
            command.append('afterok:{}'.format(':'.join(dependents_ids)))
        command.append(script_name)
        return command

    def submit_command(self, target, script_name, dependents_ids):
        self._write_script(target)

        command = self._build_submit_command(target, script_name, dependents_ids)
        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        job_id = qsub.stdout.read().strip()
        return job_id

    def build_cancel_command(self, job_ids):
        return ['scancel'] + map(str, job_ids)
