import csv
import json
import logging
import os
import os.path
import subprocess
from distutils.spawn import find_executable

from . import Backend
from ..exceptions import BackendError, NoLogFoundError
from ..utils import timer

logger = logging.getLogger(__name__)


# see squeue man page under JOB STATE CODES
JOB_STATE_CODES = {
    'BF': '?',  # BOOT_FAIL
    'CA': '?',  # CANCELLED
    'CD': '?',  # COMPLETED
    'CF': 'R',  # CONFIGURING
    'CG': 'R',  # COMPLETING
    'F': '?',   # FAILED
    'NF': '?',  # NODE_FAIL
    'PD': 'Q',  # PENDING
    'PR': '?',  # PREEMPTED
    'R': 'R',   # RUNNING
    'S': 'R',   # SUSPENDED
    'TO': '?',  # TIMEOUT
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


def _find_exe(name):
    exe = find_executable(name)
    if not exe:
        raise BackendError('Could not find executable "{}".'.format(name))
    return exe


def dump_atomic(obj, path):
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
    """

    name = 'slurm'

    supported_options = set(OPTION_TABLE.keys())
    option_defaults = {
        'cores': 1,
        'memory': '1g',
        'walltime': '01:00:00',
    }

    _JOB_DB_PATH = '.gwf/slurm-backend-jobdb.json'
    _JOB_HISTORY_PATH = '.gwf/slurm-backend-jobhistory.json'

    def configure(self, *args, **kwargs):
        self.squeue = _find_exe('squeue')
        self.sbatch = _find_exe('sbatch')
        self.scancel = _find_exe('scancel')

        # Make sure that directory for log files exists.
        self._log_dir = os.path.join(self.workflow.working_dir, '.gwf/logs')
        if not os.path.exists(self._log_dir):
            logger.debug(
                'Log directory "%s" does not exist. Creating.',
                self._log_dir
            )
            os.makedirs(self._log_dir)

        # TODO: maybe use some sort of actual db instead of a file?
        self._job_db = _read_json(self._JOB_DB_PATH)
        self._job_history = _read_json(self._JOB_HISTORY_PATH)

        self._live_job_states = self._get_live_job_states()
        logger.debug('found %d jobs', len(self._live_job_states))

        with timer('filtering jobs took %.2f ms', logger=logger):
            self._job_db = {
                target_name: self._live_job_states[job_id]
                for target_name, job_id in self._job_db.items()
                if job_id in self._live_job_states
            }

    def close(self):
        # TODO: error handling
        dump_atomic(self._job_db, self._JOB_DB_PATH)
        dump_atomic(self._job_history, self._JOB_HISTORY_PATH)

    def submitted(self, target):
        return target.name in self._job_db

    def running(self, target):
        target_job_id = self._job_db.get(target.name, None)
        return self._live_job_states.get(target_job_id, '?') == 'R'

    def submit(self, target):
        dependency_ids = [
            self._job_db[dep.name]
            for dep in self.workflow.dependencies[target]
            if dep.name in self._job_db
        ]

        cmd = [self.sbatch, "--parsable"]
        if dependency_ids:
            cmd.append(
                "--dependency=afterok:{}".format(",".join(dependency_ids))
            )

        script_contents = self._compile_script(target)
        logger.debug('Compiled script for target %s:\n%s',
                     target.name, script_contents)

        proc = subprocess.Popen(cmd)
        new_job_id, error_text = proc.communicate(script_contents)
        if proc.returncode != 0:
            raise BackendError(error_text)

        new_job_id = new_job_id.strip()

        self._job_db[target.name] = new_job_id

        if target.name not in self._job_history:
            self._job_history[target.name] = []
        self._job_history[target.name].append(new_job_id)

        # New jobs are assumed to be on-hold until the next time gwf is
        # invoked
        self._live_job_states[new_job_id] = 'H'

    def logs(self, target, stderr=False, rewind=0):
        target_history = self._job_history.get(target.name, [])
        if len(target_history) <= rewind:
            raise NoLogFoundError(
                'No log exists at this point in time.'
            )

        job_id = target_history[-1 - rewind]
        stdout_path = self._get_stdout_log_path(target, job_id=job_id)
        if stderr:
            stderr_path = self._get_stderr_log_path(target, job_id=job_id)
            return (open(stdout_path), open(stderr_path))
        return open(stdout_path)

    def cancel(self, target):
        """Cancel a target."""
        target_job_id = self._job_db.get(target.name, '?')
        if target_job_id in self._live_job_states:
            proc = subprocess.Popen([self.scancel, "-j", target_job_id])
            proc.communicate()
            # TODO: error handling

    def _get_stdout_log_path(self, target, job_id='%j'):
        return os.path.join(
            self._log_dir,
            '{}.{}.stdout'.format(target.name, job_id)
        )

    def _get_stderr_log_path(self, target, job_id='%j'):
        return os.path.join(
            self._log_dir,
            '{}.{}.stderr'.format(target.name, job_id)
        )

    def _compile_script(self, target):
        option_str = "#SBATCH {0}{1}"

        out = []
        out.append('#!/bin/bash')
        out.append('')
        out.append('####################')
        out.append('# Generated by GWF #')
        out.append('####################')
        out.append('')

        for option_name, option_value in target.options.items():
            if option_name not in self.supported_options:
                continue

            out.append(
                option_str.format(
                    OPTION_TABLE[option_name],
                    option_value
                )
            )

        out.append(option_str.format(
            '--output=', self._get_stdout_log_path(target)))
        out.append(option_str.format(
            '--error=', self._get_stderr_log_path(target)))

        out.append('')
        out.append('cd {}'.format(target.working_dir))
        out.append('export GWF_JOBID=$SLURM_JOBID')
        out.append('export GWF_TARGET_NAME="{}"'.format(target.name))
        out.append('set -e')
        out.append('')
        out.append(target.spec)
        return '\n'.join(out)

    @timer('fetching slurm job states with sqeueu took %0.2fms', logger=logger)
    def _get_live_job_states(self):
        """Ask Slurm for the state of all live jobs.

        There are two reasons why we ask for all the jobs:

            1. We don't want to spawn a subprocesses for each job
            2. Asking for multiple specific jobs would involve building a
               potentially very long commandline - which could fail if too long.

        :return: a dict mapping from job id to either 'R', 'H', 'Q' or '?'.
        """

        # The only reason this is a method instead of a function outside the
        # class is that the we need access to self.squeue. Maybe this should
        cmd = [self.squeue, '--noheader', '--format=%i;%t;%E']

        stat = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        reader = csv.reader(
            stat.stdout,
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
