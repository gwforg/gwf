"""Wrappers around code for the grid queue backend."""

import subprocess


class TorqueBackend(object):
    """Backend functionality for torque."""

    def __init__(self):
        pass

    def get_state_of_jobs(job_ids):
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


class SlurmBackend(object):
    """Backend functionality for torque."""

    def __init__(self):
        pass

    def get_state_of_jobs_slurm(job_ids):
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
