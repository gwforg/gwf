import csv
import logging
import os
import os.path
import subprocess
from distutils.spawn import find_executable

from .base import PersistableDict
from . import Backend
from ..exceptions import BackendError
from ..utils import cache, timer

logger = logging.getLogger(__name__)


# see squeue man page under JOB STATE CODES
JOB_STATE_CODES = {
    'BF': 'F',  # BOOT_FAIL
    'CA': 'F',  # CANCELLED
    'CD': 'C',  # COMPLETED
    'CF': 'R',  # CONFIGURING
    'CG': 'R',  # COMPLETING
    'F': 'F',   # FAILED
    'NF': 'F',  # NODE_FAIL
    'PD': 'Q',  # PENDING
    'PR': 'F',  # PREEMPTED
    'R': 'R',   # RUNNING
    'S': 'R',   # SUSPENDED
    'TO': 'F',  # TIMEOUT
    'SE': 'Q',  # SPECIAL_EXIT
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
    return _call_generic('squeue', '--noheader', '--format=%i;%t;%E')


def _call_scancel(job_id):
    return _call_generic('scancel_', '-j', job_id)


def _call_sbatch(script, dependencies):
    args = ['--parsable']
    if dependencies:
        args.append('--dependency=afterok:{}'.format(':'.join(dependencies)))
    return _call_generic('sbatch', *args, input=script)


@timer('Fetched job queue in %0.2fms', logger=logger)
def _get_status():
    """Ask Slurm for the state of all live jobs.

    There are two reasons why we ask for all the jobs:

        1. We don't want to spawn a subprocesses for each job
        2. Asking for multiple specific jobs would involve building a
           potentially very long commandline - which could fail if too long.

    :return: a dict mapping from job id to either 'R', 'H', 'Q' or '?'.
    """

    # The only reason this is a method instead of a function outside the
    # class is that the we need access to self.squeue. Maybe this should
    stdout, _ = _call_squeue()

    reader = csv.reader(
        stdout.splitlines(),
        delimiter=';',
        quoting=csv.QUOTE_NONE,
    )

    result = {}
    for job_id, state, depends in reader:
        simple_state = JOB_STATE_CODES[state]
        if simple_state == 'Q' and depends:
            result[job_id] = 'H'
        else:
            result[job_id] = simple_state
    return result


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

    supported_options = set(OPTION_TABLE.keys())
    option_defaults = {
        'cores': 1,
        'memory': '1g',
        'walltime': '01:00:00',
    }

    _JOB_DB_PATH = '.gwf/slurm-backend-jobdb.json'

    def __init__(self, working_dir):
        super().__init__(working_dir)
        self._tracked = PersistableDict(os.path.join(self.working_dir, self._JOB_DB_PATH))
        self._status = _get_status()

    def close(self):
        self._tracked.persist()

    def submitted(self, target):
        return target.name in self._tracked and self._tracked[target.name] in self._status

    def running(self, target):
        target_job_id = self._tracked.get(target.name, None)
        return self._status.get(target_job_id, '?') == 'R'

    def completed(self, target):
        target_job_id = self._tracked.get(target.name, None)
        return self._status.get(target_job_id, '?') == 'C'

    def submit(self, target, dependencies):
        dependency_ids = [
            self._tracked[dep.name]
            for dep in dependencies
            if dep.name in self._tracked
        ]

        script = self._compile_script(target)

        stdout, _ = _call_sbatch(script, dependency_ids)
        new_job_id = stdout.strip()

        self._tracked[target.name] = new_job_id

        # New jobs are assumed to be on-hold until the next time gwf is invoked.
        self._status[new_job_id] = 'H'

    def cancel(self, target):
        """Cancel a target."""
        job_id = self._tracked.get(target.name, '?')
        if job_id not in self._status:
            raise BackendError('Cannot cancel non-running target.')
        _call_scancel(job_id)

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
            out.append(
                option_str.format(
                    OPTION_TABLE[option_name],
                    option_value
                )
            )

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


class Script:

    def __init__(self):
        self._header_lines = []
        self._envvar_lines = []
        self._option_lines = []
        self._command = ''

    def add_header(self, text):
        self._header_lines.append('# {}'.format(text))

    def add_envvar(self, name, value):
        self._envvar_lines.append('export {}={}'.format(name, value))

    def add_option(self, option, value):
        pass

    def set_command(self, command):
        pass

    def render_script(self):
        return '\n'.join(self.header_lines + self._option_lines + self._envvar_lines + [self._command])
