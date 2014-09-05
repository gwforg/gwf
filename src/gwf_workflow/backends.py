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
                stat = subprocess.Popen(['qstat', '-f', job_id], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                for line in stat.stdout:
                    line = line.strip()
                    if line.startswith('job_state'):
                        status = line.split()[2]
                        result[job_id] = status
            except:
                pass
        return result

    def write_script_header(self, f, options):
        print >> f, '#PBS -l nodes={}:ppn={}'.format(options['nodes'], options['cores'])
        print >> f, '#PBS -l mem={}'.format(options['memory'])
        print >> f, '#PBS -l walltime={}'.format(options['walltime'])

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
            command.append('depend=afterok:{}'.format(':'.join(dependents_ids)))
        command.append(script_name)
        return command



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
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
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



AVAILABLE_BACKENDS = {
    'torque': TorqueBackend,
    'slurm': SlurmBackend,
}
