"""Wrappers around code for the grid queue backend."""

import subprocess
import os.path


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
        command = ['qsub', '-N', target.name, 
                   '-o', os.path.join(target.working_dir, 'gwf_log', target.name+'.stdout'),
                   '-e', os.path.join(target.working_dir, 'gwf_log', target.name+'.stderr'),
        if len(dependents_ids) > 0:
            command.append('-W')
            command.append('depend=afterok:{}'.format(':'.join(dependents_ids)))
        command.append(script_name)
        return command



class SlurmBackend(object):
    """Backend functionality for torque."""

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
                                     '--format=%i;%t',
                                     '--jobs', ",".join(job_ids)],
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            for line in stat.stdout:
                job_id, state = line.strip().split(';')
                result[job_id] = map_state[state]
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
        command = ['qsub', '-N', target.name, 
                   '-o', os.path.join(target.working_dir, 'gwf_log', target.name+'.stdout'),
                   '-e', os.path.join(target.working_dir, 'gwf_log', target.name+'.stderr'),
                   ]
        if len(dependents_ids) > 0:
            command.append('-W')
            command.append('depend=afterok:{}'.format(':'.join(dependents_ids)))
        command.append(script_name)
        return command



AVAILABLE_BACKENDS = {
    'torque': TorqueBackend,
    'slurm': SlurmBackend,
}