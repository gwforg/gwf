import csv
import json
import logging
import os
import os.path
import subprocess
from distutils.spawn import find_executable

from . import Backend, FileLogsMixin
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
}


@cache
def _find_exe(name):
    exe = find_executable(name)
    if exe is None:
        raise BackendError('Could not find executable "{}".'.format(name))
    return exe


def _dump_atomic(obj, path):
    with open(path + '.new', 'w') as fileobj:
        json.dump(obj, fileobj)
        fileobj.flush()
        os.fsync(fileobj.fileno())
        fileobj.close()
    os.rename(path + '.new', path)


def _read_json(path):
    try:
        with open(path) as fileobj:
            return json.load(fileobj)
    except (OSError, ValueError):
        # Catch ValueError for compatibility with Python 3.4.2. I haven't been
        # able to figure out what is different between 3.4.2 and 3.5 that
        # causes this. Essentially, 3.4.2 raises a ValueError saying that it
        # cannot parse the empty string instead of raising an OSError
        # (FileNotFoundError does not exist in 3.4.2) saying that the file does
        # not exist.
        return {}


def _call_squeue():
    squeue_path = _find_exe('squeue')
    cmd = [squeue_path, '--noheader', '--format=%i;%t;%E']

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise BackendError(stderr)
    return stdout, stderr


def _call_scancel(job_id):
    scancel_path = _find_exe('scancel')
    proc = subprocess.Popen(
        [scancel_path, "-j", job_id],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise BackendError(stderr)


def _call_sbatch(script, dependencies):
    sbatch_path = _find_exe('sbatch')

    cmd = [sbatch_path, "--parsable"]
    if dependencies:
        cmd.append(
            "--dependency=afterok:{}".format(",".join(dependencies))
        )

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = proc.communicate(script)
    if proc.returncode != 0:
        raise BackendError(stderr)
    return stdout, stderr


@timer('Fetching slurm job states with squeue took %0.2fms', logger=logger)
def _get_live_job_states():
    """Ask Slurm for the state of all live jobs.

    There are two reasons why we ask for all the jobs:

        1. We don't want to spawn a subprocesses for each job
        2. Asking for multiple specific jobs would involve building a
           potentially very long commandline - which could fail if too long.

    :return: a dict mapping from job id to either 'R', 'H', 'Q' or '?'.
    """

    # The only reason this is a method instead of a function outside the
    # class is that the we need access to self.squeue. Maybe this should
    stdout, stderr = _call_squeue()

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


class SlurmBackend(FileLogsMixin, Backend):
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
    """

    supported_options = set(OPTION_TABLE.keys())
    option_defaults = {
        'cores': 1,
        'memory': '1g',
        'walltime': '01:00:00',
    }

    _JOB_DB_PATH = '.gwf/slurm-backend-jobdb.json'

    def configure(self, working_dir, config):
        super().configure(working_dir, config)

        # Make sure that directory for log files exists.
        self._log_dir = os.path.join(self.working_dir, '.gwf/logs')
        if not os.path.exists(self._log_dir):
            logger.debug(
                'Log directory "%s" does not exist. Creating.',
                self._log_dir
            )
            os.makedirs(self._log_dir)

        # TODO: maybe use some sort of actual db instead of a file?
        self._job_db = _read_json(os.path.join(working_dir, self._JOB_DB_PATH))

        self._live_job_states = _get_live_job_states()
        logger.debug('found %d jobs', len(self._live_job_states))

        with timer('filtering jobs took %.2f ms', logger=logger):
            self._job_db = {
                target_name: job_id
                for target_name, job_id in self._job_db.items()
                if job_id in self._live_job_states
            }

    def close(self):
        if hasattr(self, '_job_db'):
            _dump_atomic(
                self._job_db,
                os.path.join(self.working_dir, self._JOB_DB_PATH)
            )

    def submitted(self, target):
        return target.name in self._job_db

    def running(self, target):
        target_job_id = self._job_db.get(target.name, None)
        return self._live_job_states.get(target_job_id, '?') == 'R'

    def failed(self, target):
        target_job_id = self._job_db.get(target.name, None)
        return self._live_job_states.get(target_job_id, '?') == 'F'

    def completed(self, target):
        target_job_id = self._job_db.get(target.name, None)
        return self._live_job_states.get(target_job_id, '?') == 'C'

    def submit(self, target, dependencies):
        dependency_ids = [
            self._job_db[dep.name]
            for dep in dependencies
            if dep.name in self._job_db
        ]

        script = self._compile_script(target)

        stdout, stderr = _call_sbatch(script, dependency_ids)
        new_job_id = stdout.strip()

        self._job_db[target.name] = new_job_id

        # New jobs are assumed to be on-hold until the next time gwf is
        # invoked
        self._live_job_states[new_job_id] = 'H'

    def cancel(self, target):
        """Cancel a target."""
        job_id = self._job_db.get(target.name, '?')
        if job_id not in self._live_job_states:
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

        out.append(option_str.format(
            '--output=', FileLogsMixin.stdout_path(target)))
        out.append(option_str.format(
            '--error=', FileLogsMixin.stderr_path(target)))

        out.append('')
        out.append('cd {}'.format(target.working_dir))
        out.append('export GWF_JOBID=$SLURM_JOBID')
        out.append('export GWF_TARGET_NAME="{}"'.format(target.name))
        out.append('set -e')
        out.append('')
        out.append(target.spec)
        return '\n'.join(out)
