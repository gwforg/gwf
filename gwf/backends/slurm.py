"""This modules implements a GWF backend for the Slurm workload manager."""

from __future__ import absolute_import, print_function

import distutils.spawn
import json
import logging
import subprocess

from ..core import PreparedWorkflow
from ..exceptions import BackendError
from ..utils import cache, dfs, timer
from .base import Backend

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


def _find_exe(name):
    return distutils.spawn.find_executable(name)


SQUEUE = _find_exe('squeue')
SBATCH = _find_exe('sbatch')
SCANCEL = _find_exe('scancel')


@cache
@timer('fetching slurm job states with sqeueu took %0.2fms', logger=logger)
def _live_job_states():
    """Ask Slurm for the state of all live jobs.

    There are two reasons why we ask for all the jobs:

        1. We don't want to spawn a subprocesses for each job
        2. Asking for multiple specific jobs would involve building a
           potentially very long commandline - which could fail if too long.

    :return: a dict mapping from job id to either R, H or Q.
    """
    cmd = [SQUEUE, '--noheader', '--format=%i;%t;%E']
    stat = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)

    result = {}
    for line in stat.stdout:
        job_id, state, depends = line.strip().split(';')
        simple_state = JOB_STATE_CODES[state]
        if simple_state == 'Q' and depends != '':
            result[job_id] = 'H'
        else:
            result[job_id] = simple_state
    return result


def _compile_script(target):
    out = []
    out.append('#!/bin/bash')

    option_table = [
        ("-N ", "nodes"),
        ("-c ", "cores"),
        ("--mem=", "memory"),
        ("-t ", "walltime"),
        ("-p ", "queue"),
        ("-A ", "account"),
        ("-C ", "constraint"),
        ("--mail-type=", "mail_type"),
        ("--mail-user=", "mail_user"),
    ]

    out.extend(
        "#SBATCH {0}{1}".format(slurm_flag, target.options[gwf_name])
        for slurm_flag, gwf_name in option_table
        if gwf_name in target.options
    )

    out.append('')
    out.append('# GWF generated code ...')
    out.append('cd {}'.format(target.working_dir))
    out.append('export GWF_JOBID=$SLURM_JOBID')
    out.append('set -e')
    out.append('')

    out.append('# Script from workflow')
    out.append(target.spec)
    return '\n'.join(out)


class SlurmBackend(Backend):
    """Backend for the slurm workload manager."""

    name = 'slurm'

    def __init__(self, workflow):
        super().__init__(workflow)

        # TODO: maybe use some sort of actual db instead of a file?
        try:
            with open(".gwf/slurm-backend-jobdb.json") as fileobj:
                self.job_db = json.load(fileobj)
        except FileNotFoundError:
            self.job_db = {}

        self.live_job_states = _live_job_states()
        logger.debug('found %d jobs', len(self.live_job_states))

        with timer('filtering jobs took %.2f ms', logger=logger):
            for target_name, job_id in self.job_db.items():
                if job_id not in self.live_job_states:
                    del self.job_db[target_name]

    def close(self):
        # TODO: error handling
        with open(".gwf/slurm-backend-jobdb.json", "w") as fileobj:
            json.dump(self.job_db, fileobj)

    def submitted(self, target):
        """Return whether the target has been submitted."""
        return target.name in self.job_db

    def running(self, target):
        """Return whether the target is running."""
        target_job_id = self.job_db.get(target.name, None)
        return self.live_job_states.get(target_job_id, '?') == 'R'

    def submit(self, target):
        """Submit a target."""
        dependency_ids = [
            self.job_db[dep.name]
            for dep in self.workflow.dependencies[target]
            if dep.name in self.job_db
        ]

        cmd = [SBATCH, "--parsable"]
        if dependency_ids:
            cmd.append(
                "--dependency=afterok:{}".format(",".join(dependency_ids))
            )

        script_contents = _compile_script(target)

        proc = subprocess.Popen(cmd)
        new_job_id, error_text = proc.communicate(script_contents)
        new_job_id = new_job_id.strip()
        if proc.returncode != 0:
            raise BackendError(error_text)

        self.job_db[target.name] = new_job_id
        # New jobs are assumed to be on-hold until the next time gwf is invoked
        self.live_job_states[new_job_id] = 'H'

    def cancel(self, target):
        """Cancel a target."""
        target_job_id = self.job_db.get(target.name, '?')
        if target_job_id in self.live_job_states:
            proc = subprocess.Popen([SCANCEL, "-j", target_job_id])
            stdout, stderr = proc.communicate()
            # TODO: error handling
