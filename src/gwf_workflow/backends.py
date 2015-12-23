"""Wrappers around code for the grid queue backend."""

import subprocess
import os
import os.path


def _mkdir_if_not_exist(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


class TorqueBackend(object):
    """Backend functionality for torque."""

    def __init__(self):
        pass

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
        print >> f, '#PBS -l nodes={}:ppn={}'.format(
            options['nodes'], options['cores'])
        print >> f, '#PBS -l mem={}'.format(options['memory'])
        print >> f, '#PBS -l walltime={}'.format(options['walltime'])
        if 'queue' in options:
            print >> f, '#PBS -q {}'.format(options['queue'])

    def write_script_variables(self, f):
        print >> f, 'export GWF_JOBID=$PBS_JOBID'

    def build_submit_command(self, target, script_name, dependents_ids):
        log_dir = os.path.join(target.working_dir, 'gwf_log')
        _mkdir_if_not_exist(log_dir)

        command = ['qsub', '-N', target.name,
                   '-o', os.path.join(log_dir, target.name+'.stdout'),
                   '-e', os.path.join(log_dir, target.name+'.stderr'),
                   ]
        if len(dependents_ids) > 0:
            command.append('-W')

            ids = ':'.join(dependents_ids)
            command.append('depend=afterok:{}'.format(ids))
        command.append(script_name)
        return command

    def submit_command(self, target, script_name, dependents_ids):
        command = self.build_submit_command(target, script_name, dependents_ids)
        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        job_id = qsub.stdout.read().strip()
        return job_id

    def build_cancel_command(self, job_ids):
        return ['qdel'] + map(str, job_ids)

class SlurmBackend(object):
    """Backend functionality for slurm."""

    def __init__(self):
        pass

    def get_state_of_jobs(self, job_ids):
        result = dict((job_id, False) for job_id in job_ids)
        map_state = {  # see squeue man page under JOB STATE CODES
                       'BF': '?',  # BOOT_FAIL
                       'CA': '?',  # CANCELLED
                       'CD': '?',  # COMPLETED
                       'CF': 'R',  # CONFIGURING
                       'CG': 'R',  # COMPLETING
                       'F':  '?',  # FAILED
                       'NF': '?',  # NODE_FAIL
                       'PD': 'Q',  # PENDING
                       'PR': '?',  # PREEMPTED
                       'R':  'R',  # RUNNING
                       'S':  'R',  # SUSPENDED
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
        print >> f, '#SBATCH -N {}'.format(options['nodes'])
        print >> f, '#SBATCH -c {}'.format(options['cores'])
        print >> f, '#SBATCH --mem={}'.format(options['memory'])
        print >> f, '#SBATCH -t {}'.format(options['walltime'])
        if 'queue' in options:
            print >> f, '#SBATCH -p {}'.format(options['queue'])

    def write_script_variables(self, f):
        print >> f, 'export GWF_JOBID=$SLURM_JOBID'

    def build_submit_command(self, target, script_name, dependents_ids):
        log_dir = os.path.join(target.working_dir, 'gwf_log')
        _mkdir_if_not_exist(log_dir)

        command = ['sbatch', '-J', target.name, '--parsable',
                   '-o', os.path.join(log_dir, target.name+'.%j.stdout'),
                   '-e', os.path.join(log_dir, target.name+'.%j.stderr'),
                   ]
        if len(dependents_ids) > 0:
            command.append('-d')
            command.append('afterok:{}'.format(':'.join(dependents_ids)))
        command.append(script_name)
        return command

    def submit_command(self, target, script_name, dependents_ids):
        command = self.build_submit_command(target, script_name, dependents_ids)
        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        job_id = qsub.stdout.read().strip()
        return job_id

    def build_cancel_command(self, job_ids):
        return ['scancel'] + map(str, job_ids)

class LocalBackend(object):
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

    def build_submit_command(self, target, script_name, dependent_ids):
        command = ["bash", script_name]
        return command

    def submit_command(self, target, script_name, dependents_ids):
        log_dir = os.path.join(target.working_dir, 'gwf_log')
        command = self.build_submit_command(target, script_name, dependents_ids)
        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(os.path.join(log_dir, target.name+'.stdout'), "w") as outfile:
            print >> outfile, qsub.stdout.read()
        with open(os.path.join(log_dir, target.name+'.stderr'), "w") as outfile:
            print >> outfile, qsub.stderr.read()
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
