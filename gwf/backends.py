from __future__ import absolute_import
from __future__ import print_function

import subprocess
import os
import os.path

import gwf

from gwf.helpers import make_absolute_path, escape_file_name


# This will be set in the gwf script and refer to the grid backend used.
BACKEND = None


def _mkdir_if_not_exist(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


class Backend(object):

    def __init__(self, working_dir):
        pass

    def get_state_of_jobs(self):
        raise NotImplementedError()

    def write_script_header(self, f, options):
        raise NotImplementedError()

    def write_script_variables(self, f):
        raise NotImplementedError()

    def submit_command(self, target, script_name, dependent_ids):
        raise NotImplementedError()

    def build_cancel_command(self, job_ids):
        raise NotImplementedError()


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

        log_dir = os.path.join(target.working_dir, 'gwf_log')
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


class SlurmBackend(Backend):
    """Backend functionality for slurm."""

    def __init__(self):

        self.script_dir = make_absolute_path(gwf.WORKING_DIR, '.scripts')
        self.script_name = make_absolute_path(self.script_dir, escape_file_name(self.name))

    def make_script_dir(self):
        script_dir = self.script_dir
        # Don't use the _file_exists() function here. It caches its status and that won't work for the script dir.
        if not os.path.exists(script_dir):
            os.makedirs(script_dir)

    def write_script(self, target):
        """Write the code to a script that can be executed."""

        self.make_script_dir()
        f = open(self.script_name, 'w')

        print("#!/bin/bash", file=f)

        self.write_script_header(f, target.options)
        print(file=f)

        print('# GWF generated code ...', file=f)
        print('cd %s' % gwf.WORKING_DIR, file=f)
        self.write_script_variables(f)
        print(f, "set -e", file=f)
        print(file=f)

        print('# Script from workflow', file=f)

        print(self.spec, file=f)

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

    def write_script_header(self, f, options):
        print('#SBATCH -N {}'.format(options['nodes']), file=f)
        print('#SBATCH -c {}'.format(options['cores']), file=f)
        print('#SBATCH --mem={}'.format(options['memory']), file=f)
        print('#SBATCH -t {}'.format(options['walltime']), file=f)
        if 'queue' in options:
            print('#SBATCH -p {}'.format(options['queue']), file=f)
        if 'account' in options:
            print('#SBATCH -A {}'.format(options['account']), file=f)
        if 'constraint' in options:
            print('#SBATCH -C {}'.format(options['constraint']), file=f)
        if 'mail_type' in options:
            print('#SBATCH --mail-type={}'.format(options['mail_type']), file=f)
        if 'mail_user' in options:
            print('#SBATCH --mail-user={}'.format(options['mail_user']), file=f)

    def write_script_variables(self, f):
        print('export GWF_JOBID=$SLURM_JOBID', file=f)

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
        self.write_script(target)

        command = self._build_submit_command(target, script_name, dependents_ids)
        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        job_id = qsub.stdout.read().strip()
        return job_id

    def build_cancel_command(self, job_ids):
        return ['scancel'] + map(str, job_ids)


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
        print >> f, 'export GWF_JOBID={}'.format(self.next_job_id)

    def _build_submit_command(self, target, script_name, dependent_ids):
        command = ["bash", script_name]
        return command

    def submit_command(self, target, script_name, dependents_ids):
        log_dir = os.path.join(target.working_dir, 'gwf_log')
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


AVAILABLE_BACKENDS = {
    'torque': TorqueBackend,
    'slurm': SlurmBackend,
    'local': LocalBackend
}
