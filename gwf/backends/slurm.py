import csv
import logging
import subprocess
from distutils.spawn import find_executable

from .base import PersistableDict, Status, UnknownDependencyError, UnknownTargetError
from . import Backend
from ..exceptions import BackendError
from ..utils import cache, timer

logger = logging.getLogger(__name__)


UNKNOWN_DEPENDENCY_EXC = 'Cannot submit {} with dependency {} since this target has not been submitted yet.'


JOB_STATE_CODES = {
    'BF': Status.UNKNOWN,  # BOOT_FAIL
    'CA': Status.UNKNOWN,  # CANCELLED
    'CD': Status.UNKNOWN,  # COMPLETED
    'CF': Status.RUNNING,  # CONFIGURING
    'CG': Status.RUNNING,  # COMPLETING
    'F': Status.UNKNOWN,   # FAILED
    'NF': Status.UNKNOWN,  # NODE_FAIL
    'PD': Status.SUBMITTED,  # PENDING
    'PR': Status.UNKNOWN,  # PREEMPTED
    'R': Status.RUNNING,   # RUNNING
    'S': Status.RUNNING,   # SUSPENDED
    'TO': Status.UNKNOWN,  # TIMEOUT
    'SE': Status.SUBMITTED,  # SPECIAL_EXIT
}


OPTION_TABLE = {
    "nodes": "-N ",
    "cores": "-c ",
    "memory": "--mem=",
    "walltime": "-t ",
    "queue": "-p ",
    "account": "-A ",
    "constraint": "-C ",
    "mail_type": "--mail-type=",
    "mail_user": "--mail-user=",
    "qos": "--qos="
}


@cache
def _find_exe(name):
    exe = find_executable(name)
    if exe is None:
        raise BackendError('Could not find executable "{}".'.format(name))
    return exe


def _call_generic(executable_name, *args, input=None):
    executable_path = _find_exe(executable_name)
    proc = subprocess.Popen(
        [executable_path] + list(args),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = proc.communicate(input)
    if proc.returncode != 0:
        raise BackendError(stderr)
    return stdout, stderr


def _call_sacct(job_id):
    return _call_generic('sacct', '--noheader', '--long', '--parsable2', '--allocations', '--jobs', job_id)


def _call_squeue():
    return _call_generic('squeue', '--noheader', '--format=%i;%t')


def _call_scancel(job_id):
    return _call_generic('scancel', '-j', job_id)


def _call_sbatch(script, dependencies):
    args = ['--parsable']
    if dependencies:
        args.append('--dependency=afterok:{}'.format(':'.join(dependencies)))
    return _call_generic('sbatch', *args, input=script)


@timer('Fetched job queue in %0.2fms.', logger=logger)
def init_status_from_queue():
    """Ask Slurm for the state of all live jobs.

    There are two reasons why we ask for all the jobs:

        1. We don't want to spawn a subprocesses for each job
        2. Asking for multiple specific jobs would involve building a
           potentially very long commandline - which could fail if too long.
    """
    stdout, _ = _call_squeue()
    job_states = {}
    for line in stdout.splitlines():
        job_id, state = line.split(';')
        job_states[job_id] = JOB_STATE_CODES[state]
    return job_states


class SlurmBackend(Backend):
    """Backend for the Slurm workload manager.

    To use this backend you must activate the `slurm` backend.

    **Backend options:**

    None available.

    **Target options:**

    * **cores (int):**
      Number of cores allocated to this target (default: 1).
    * **memory (str):**
      Memory allocated to this target (default: 1).
    * **walltime (str):**
      Time limit for this target (default: 01:00:00).
    * **queue (str):**
      Queue to submit the target to. To specify multiple queues, specify a
      comma-separated list of queue names.
    * **account (str):**
      Account to be used when running the target.
    * **constraint (str):**
      Constraint string. Equivalent to setting the `--constraint` flag on
      `sbatch`.
    * **qos (str):**
      Quality-of-service strring. Equivalent to setting the `--qos` flog
      on `sbatch`.
    """

    option_defaults = {
        'cores': 1,
        'memory': '1g',
        'walltime': '01:00:00',
        'nodes': None,
        'queue': None,
        'account': None,
        'constraint': None,
        'mail_type': None,
        'mail_user': None,
        'qos': None,
    }

    def __init__(self):
        self._tracked = PersistableDict('.gwf/slurm-backend-tracked.json')
        self._status = init_status_from_queue()

        for job_name, job_id in list(self._tracked.items()):
            if job_id not in self._status:
                del self._tracked[job_name]

    def status(self, target):
        try:
            return self.get_status(target)
        except KeyError:
            return Status.UNKNOWN

    def submit(self, target, dependencies):
        script = self._compile_script(target)
        dependency_ids = self.collect_dependency_ids(dependencies)
        stdout, _ = _call_sbatch(script, dependency_ids)
        job_id = stdout.strip()
        self.add_job(target, job_id)

    def cancel(self, target):
        """Cancel a target."""
        try:
            job_id = self.get_job_id(target)
        except KeyError:
            raise UnknownTargetError(target.name)

        _call_scancel(job_id)
        self.forget_job(target)

    def close(self):
        self._tracked.persist()

    def _compile_script(self, target):
        option_str = "#SBATCH {0}{1}"

        out = []
        out.append('#!/bin/bash')
        out.append('')
        out.append('####################')
        out.append('# Generated by GWF #')
        out.append('####################')
        out.append('')

        out.append(option_str.format('--job-name=', target.name))

        for option_name, option_value in target.options.items():
            out.append(option_str.format(OPTION_TABLE[option_name], option_value))

        out.append(option_str.format('--output=', self._log_path(target, 'stdout')))
        out.append(option_str.format('--error=', self._log_path(target, 'stderr')))

        out.append('')
        out.append('cd {}'.format(target.working_dir))
        out.append('export GWF_JOBID=$SLURM_JOBID')
        out.append('export GWF_TARGET_NAME="{}"'.format(target.name))
        out.append('set -e')
        out.append('')
        out.append(target.spec)
        return '\n'.join(out)

    def add_job(self, target, job_id, initial_status=Status.SUBMITTED):
        self.set_job_id(target, job_id)
        self.set_status(target, initial_status)

    def forget_job(self, target):
        job_id = self.get_job_id(target)
        del self._status[job_id]
        del self._tracked[target.name]

    def get_job_id(self, target):
        return self._tracked[target.name]

    def set_job_id(self, target, job_id):
        self._tracked[target.name] = job_id

    def get_status(self, target):
        job_id = self.get_job_id(target)
        return self._status[job_id]

    def set_status(self, target, status):
        job_id = self.get_job_id(target)
        self._status[job_id] = status

    def collect_dependency_ids(self, dependencies):
        try:
            return [self._tracked[dep.name] for dep in dependencies]
        except KeyError as exc:
            raise UnknownDependencyError(exc.args[0])
