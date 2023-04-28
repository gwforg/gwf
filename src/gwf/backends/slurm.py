"""Backend for the Slurm workload manager.

To use this backend you must activate the `slurm` backend.

**Backend options:**

* **backend.slurm.log_mode (str):** Must be either `full`, `merged` or
    `none`. If `full`, two log files will be stored for each target, one for
    standard output and one for standard error. If `merged`, only one log
    file will be written containing the combined streams. If `none`, no logs
    will be stored. (default: `full`).

* **backend.slurm.accounting_enabled (str):** If enabled, will use `sacct` to
  fetch job status from Slurm. This enables *gwf* to report when targets have
  failed (default: true).

**Target options:**

* **cores (int):**
    Number of cores allocated to this target (default: 1).
* **memory (str):**
    Memory allocated to this target (default: 1).
* **walltime (str):**
    Time limit for this target (default: 01:00:00).
* **queue (str):**
    Queue to submit the target to. To specify multiple queues, specify a
    comma-separated list of queue names. A queue is equivalent to a Slurm
    partition.
* **account (str):**
    Account to be used when running the target.
* **constraint (str):**
    Constraint string. Equivalent to setting the `--constraint` flag on
    `sbatch`.
* **qos (str):**
    Quality-of-service string. Equivalent to setting the `--qos` flog
    on `sbatch`.
* **mail_type (str):**
    Equivalent to the `--mail-type` flag on `sbatch`.
* **mail_user (str):**
    Equivalent to the `--mail-user` flag on `sbatch`.
* **gres (str):**
    Equivalent to the `--gres` flog on `sbatch`. Usually used to
    request access to GPUs.
"""

import logging
import os.path
from collections import defaultdict

import attrs

from ..utils import ensure_trailing_newline
from .base import BackendStatus, TrackingBackend
from .utils import call, has_exe

logger = logging.getLogger(__name__)


SLURM_SHORT_STATES = {
    "BF": BackendStatus.FAILED,  # Job terminated due to launch failure, typically due to a hardware failure (e.g. unable to boot the node or block and the job can not be requeued).
    "CA": BackendStatus.FAILED,  # Job was explicitly cancelled by the user or system administrator. The job may or may not have been initiated.
    "CD": BackendStatus.COMPLETED,  # Job has terminated all processes on all nodes with an exit code of zero.
    "CF": BackendStatus.SUBMITTED,  # Job has been allocated resources, but are waiting for them to become ready for use (e.g. booting).
    "CG": BackendStatus.RUNNING,  # Job is in the process of completing. Some processes on some nodes may still be active.
    "DL": BackendStatus.FAILED,  # Job terminated on deadline.
    "F": BackendStatus.FAILED,  # Job terminated with non-zero exit code or other failure condition.
    "NF": BackendStatus.FAILED,  # Job terminated due to failure of one or more allocated nodes.
    "OOM": BackendStatus.FAILED,  # Job experienced out of memory error.
    "PD": BackendStatus.SUBMITTED,  # Job is awaiting resource allocation.
    "PR": BackendStatus.FAILED,  # Job terminated due to preemption.
    "R": BackendStatus.RUNNING,  # Job currently has an allocation.
    "RD": BackendStatus.SUBMITTED,  # Job is being held after requested reservation was deleted.
    "RF": BackendStatus.SUBMITTED,  # Job is being requeued by a federation.
    "RH": BackendStatus.SUBMITTED,  # Held job is being requeued.
    "RQ": BackendStatus.SUBMITTED,  # Completing job is being requeued.
    "RS": BackendStatus.SUBMITTED,  # Job is about to change size.
    "RV": BackendStatus.SUBMITTED,  # Sibling was removed from cluster due to other cluster starting the job.
    "SE": BackendStatus.SUBMITTED,  # The job was requeued in a special state. This state can be set by users, typically in EpilogSlurmctld, if the job has terminated with a particular exit value.
    "SO": BackendStatus.SUBMITTED,  # Job is staging out files.
    "ST": BackendStatus.RUNNING,  # Job has an allocation, but execution has been stopped with SIGSTOP signal. CPUS have been retained by this job.
    "S": BackendStatus.RUNNING,  # Job has an allocation, but execution has been suspended and CPUs have been released for other jobs.
    "TO": BackendStatus.FAILED,  # Job terminated upon reaching its time limit.
}

SLURM_LONG_STATES = {
    "BOOT_FAIL": "BF",  # Job terminated due to launch failure, typically due to a hardware failure (e.g. unable to boot the node or block and the job can not be requeued).
    "CANCELLED": "CA",  # Job was explicitly cancelled by the user or system administrator. The job may or may not have been initiated.
    "COMPLETED": "CD",  # Job has terminated all processes on all nodes with an exit code of zero.
    "DEADLINE": "DL",  # Job terminated on deadline.
    "FAILED": "F",  # Job terminated with non-zero exit code or other failure condition.
    "NODE_FAIL": "NF",  # Job terminated due to failure of one or more allocated nodes.
    "OUT_OF_MEMORY": "OOM",  # Job experienced out of memory error.
    "PENDING": "PD",  # Job is awaiting resource allocation.
    "PREEMPTED": "PR",  # Job terminated due to preemption.
    "RUNNING": "R",  # Job currently has an allocation.
    "REQUEUED": "RQ",  # Job was requeued.
    "RESIZING": "RS",  # Job is about to change size.
    "REVOKED": "RV",  # Sibling was removed from cluster due to other cluster starting the job.
    "SUSPENDED": "S",  # Job has an allocation, but execution has been suspended and CPUs have been released for other jobs.
    "TIMEOUT": "TO",  # Job terminated upon reaching its time limit.
}


SLURM_JOB_STATES = defaultdict(lambda: BackendStatus.UNKNOWN, SLURM_SHORT_STATES)


TARGET_DEFAULTS = {
    "cores": 1,
    "memory": "1g",
    "walltime": "01:00:00",
    "nodes": None,
    "queue": None,
    "account": None,
    "constraint": None,
    "mail_type": None,
    "mail_user": None,
    "qos": None,
    "gres": None,
}

OPTION_FLAGS = {
    "nodes": "-N ",
    "cores": "-c ",
    "memory": "--mem=",
    "walltime": "-t ",
    "queue": "-p ",
    "account": "-A ",
    "constraint": "-C ",
    "mail_type": "--mail-type=",
    "mail_user": "--mail-user=",
    "qos": "--qos=",
    "gres": "--gres=",
}

OPTION_STR = "#SBATCH {0}{1}"


@attrs.define
class SlurmOps:
    working_dir: str = attrs.field()
    log_mode: str = attrs.field()
    accounting_enabled: bool = attrs.field()
    target_defaults: dict = attrs.field()

    def cancel_job(self, job_id):
        # The --verbose flag here is necessary, otherwise we're not able to tell
        # whether the command failed. See the comment in call() if you
        # want to know more.
        call("scancel", "--verbose", job_id)

    def submit_target(self, target, dependencies):
        script = self.compile_script(target)
        args = ["--parsable"]
        if dependencies:
            args.append("--dependency=afterok:{}".format(":".join(dependencies)))
        return call("sbatch", *args, input=script).strip()

    def get_job_states_from_squeue(self, tracked_jobs):
        logger.debug("Loading job states from squeue")
        job_states = {}
        for line in call(
            "squeue", "--noheader", "--format=%i;%t", "--all"
        ).splitlines():
            job_id, state = line.strip().split(";")
            job_states[job_id] = SLURM_JOB_STATES[state]
        return job_states

    def get_job_states_from_sacct(self, tracked_jobs):
        if not tracked_jobs:
            return {}

        cmd = [
            "sacct",
            "--noheader",
            "--parsable2",
            "--format=jobid,state",
            "--allocations",
            "--jobs",
            ",".join(tracked_jobs),
        ]
        job_states = {}
        for line in call(*cmd).splitlines():
            job_id, state = line.strip().split("|")
            state = SLURM_LONG_STATES[state]
            job_states[job_id] = SLURM_JOB_STATES[state]
        return job_states

    def get_job_states_from_sacct_batched(self, tracked_jobs, batch_size=1024):
        logger.debug("Loading job states from sacct")
        job_states = {}
        for idx in range(0, len(tracked_jobs), batch_size):
            batch = tracked_jobs[idx : idx + batch_size]
            job_states.update(self.get_job_states_from_sacct(batch))
        return job_states

    def get_job_states(self, tracked_jobs):
        if self.accounting_enabled:
            return self.get_job_states_from_sacct_batched(tracked_jobs)
        return self.get_job_states_from_squeue(tracked_jobs)

    def compile_script(self, target):
        out = []
        out.append("#!/bin/bash")
        out.append("# Generated by: gwf")

        out.append(OPTION_STR.format("--job-name=", target.name))

        for option_name, option_value in target.options.items():
            out.append(OPTION_STR.format(OPTION_FLAGS[option_name], option_value))

        if self.log_mode == "full":
            out.append(
                OPTION_STR.format(
                    "--output=",
                    os.path.join(
                        self.working_dir, ".gwf", "logs", target.name + ".stdout"
                    ),
                )
            )
            out.append(
                OPTION_STR.format(
                    "--error=",
                    os.path.join(
                        self.working_dir, ".gwf", "logs", target.name + ".stderr"
                    ),
                )
            )
        elif self.log_mode == "merged":
            out.append(
                OPTION_STR.format(
                    "--output=",
                    os.path.join(
                        self.working_dir, ".gwf", "logs", target.name + ".stdout"
                    ),
                )
            )
        elif self.log_mode == "none":
            out.append(OPTION_STR.format("--output=", "/dev/null"))

        out.append("")
        out.append("cd {}".format(target.working_dir))
        out.append("export GWF_JOBID=$SLURM_JOBID")
        out.append('export GWF_TARGET_NAME="{}"'.format(target.name))
        out.append("set -e")
        out.append("")
        out.append(ensure_trailing_newline(target.spec))
        return "\n".join(out)

    def close(self):
        pass


def create_backend(working_dir, log_mode="full", accounting_enabled=True):
    return TrackingBackend(
        working_dir,
        name="slurm",
        ops=SlurmOps(
            working_dir, log_mode, accounting_enabled, target_defaults=TARGET_DEFAULTS
        ),
    )


def priority():
    if has_exe("sbatch") and has_exe("sinfo"):
        return 100
    return -100


setup = (create_backend, priority())
