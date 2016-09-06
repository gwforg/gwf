"""This modules implements a GWF backend for the Slurm workload manager."""

from __future__ import absolute_import, print_function

import distutils.spawn
import json
import logging
import os
import os.path
import subprocess
import sys

from ..core import PreparedWorkflow
from ..exceptions import WorkflowNotPreparedError
from ..utils import cache, dfs
from .base import Backend

logger = logging.getLogger(__name__)


@cache
def _find_slurm_executable(name):
    return distutils.spawn.find_executable(name)


@cache
def _live_job_states():
    """Ask Slurm for the state of ALL live jobs.

    There are two reasons why we ask for all the jobs:
        1) We don't want to spawn a subprocesses for each job
        2) Asking for multiple specific jobs would involve building a
        potentially very long commandline - which could fail if too long.

    The result is a dict mapping from jobid to one of 'RQH'.
    """
    result = dict()
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
        squeue = _find_slurm_executable("squeue")
        stat = subprocess.Popen([squeue,
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


class SlurmBackend(Backend):
    """Backend for the slurm workload manager."""

    name = 'slurm'

    def __init__(self, workflow):
        super().__init__(workflow)

        # TODO: maybe use some sort of actual db instead of a file?
        if os.path.exists(".gwf/slurm-backend-jobdb.json"):
            self.job_db = json.loads(
                open(".gwf/slurm-backend-jobdb.json").read())
        else:
            self.job_db = dict()

        self.live_job_states = _live_job_states()
        for target_name, job_id in list(self.job_db.items()):
            if job_id not in self.live_job_states:
                del self.job_db[target_name]

    def close(self):
        """Close the job database"""
        encoded = json.dumps(self.job_db)
        # TODO: error handling
        with open(".gwf/slurm-backend-jobdb-new.json", "w") as f:
            f.write(encoded)
            f.write("\n")
        os.rename(".gwf/slurm-backend-jobdb-new.json",
                  ".gwf/slurm-backend-jobdb.json")

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

        sbatch = _find_slurm_executable("sbatch")
        cmd = [sbatch, "--parsable"]
        if dependency_ids:
            cmd.append(
                "--dependency=afterok:{}".format(",".join(dependency_ids))
            )
        proc = subprocess.Popen(cmd)
        script_contents = _compile_script(target)
        new_job_id, error_text = proc.communicate(script_contents)
        new_job_id = new_job_id.strip()
        if proc.returncode != 0:
            # TODO: better error handling
            print(error_text, file=sys.stderr)
            sys.exit(1)
        self.job_db[target.name] = new_job_id
        # New jobs are assumed to be on-hold until the next time gwf is invoked
        self.live_job_states[new_job_id] = 'H'

    def cancel(self, target):
        """Cancel a target."""
        target_job_id = self.job_db.get(target.name, '?')
        if target_job_id in self.live_job_states:
            scancel = _find_slurm_executable("scancel")
            proc = subprocess.Popen([scancel, "-j", target_job_id])
            stdout, stderr = proc.communicate()
            # TODO: error handling


def _compile_script(target):
    options = target.options
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
        "#SBATCH {0}{1}".format(slurm_flag, options[gwf_name])
        for slurm_flag, gwf_name in option_table
        if gwf_name in options
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
