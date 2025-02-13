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
import sys
import io
from collections import defaultdict
from contextlib import contextmanager

import attrs

from ..executors import serialize
from .base import BackendStatus, TrackingBackend
from .utils import call, has_exe

logger = logging.getLogger(__name__)


SLURM_SHORT_STATES = {
    # Job terminated due to launch failure, typically due to a hardware failure (e.g.
    # unable to boot the node or block and the job can not be requeued).
    "BF": BackendStatus.FAILED,
    # Job was explicitly cancelled by the user or system administrator. The job may or
    # may not have been initiated.
    "CA": BackendStatus.CANCELLED,
    # Job has terminated all processes on all nodes with an exit code of zero.
    "CD": BackendStatus.COMPLETED,
    # Job has been allocated resources, but are waiting for them to become ready for use
    # (e.g. booting).
    "CF": BackendStatus.SUBMITTED,
    # Job is in the process of completing. Some processes on some nodes may still be
    # active.
    "CG": BackendStatus.RUNNING,
    # Job terminated on deadline.
    "DL": BackendStatus.FAILED,
    # Job terminated with non-zero exit code or other failure condition.
    "F": BackendStatus.FAILED,
    # Job terminated due to failure of one or more allocated nodes.
    "NF": BackendStatus.FAILED,
    # Job experienced out of memory error.
    "OOM": BackendStatus.FAILED,
    # Job is awaiting resource allocation.
    "PD": BackendStatus.SUBMITTED,
    # Job terminated due to preemption.D,
    "PR": BackendStatus.FAILED,
    # Job currently has an allocation.
    "R": BackendStatus.RUNNING,
    # Job is being held after requested reservation was deleted.
    "RD": BackendStatus.SUBMITTED,
    # Job is being requeued by a federation.
    "RF": BackendStatus.SUBMITTED,
    # Held job is being requeued.
    "RH": BackendStatus.SUBMITTED,
    # Completing job is being requeued.
    "RQ": BackendStatus.SUBMITTED,
    # Job is about to change size.
    "RS": BackendStatus.SUBMITTED,
    # Sibling was removed from cluster due to other cluster starting the job.
    "RV": BackendStatus.SUBMITTED,
    # The job was requeued in a special state. This state can be set by users, typically
    # in EpilogSlurmctld, if the job has terminated with a particular exit value.
    "SE": BackendStatus.SUBMITTED,
    # Job is staging out files.
    "SO": BackendStatus.SUBMITTED,
    # Job has an allocation, but execution has been stopped with SIGSTOP signal. CPUS
    # have been retained by this job.,
    "ST": BackendStatus.RUNNING,
    # Job has an allocation, but execution has been suspended and CPUs have been
    # released for other jobs.
    "S": BackendStatus.RUNNING,
    # Job terminated upon reaching its time limit.
    "TO": BackendStatus.FAILED,
}

SLURM_LONG_STATES = {
    # Job terminated due to launch failure, typically due to a hardware failure (e.g.
    # unable to boot the node or block and the job can not be requeued).
    "BOOT_FAIL": "BF",
    # Job was explicitly cancelled by the user or system administrator. The job may or
    # may not have been initiated.
    "CANCELLED": "CA",
    # Job has terminated all processes on all nodes with an exit code of zero.
    "COMPLETED": "CD",
    # Job terminated on deadline.
    "DEADLINE": "DL",
    # Job terminated with non-zero exit code or other failure condition.,
    "FAILED": "F",
    # Job terminated due to failure of one or more allocated nodes.
    "NODE_FAIL": "NF",
    # Job experienced out of memory error.
    "OUT_OF_MEMORY": "OOM",
    # Job is awaiting resource allocation.OOM",
    "PENDING": "PD",
    # Job terminated due to preemption.
    "PREEMPTED": "PR",
    # Job currently has an allocation.,
    "RUNNING": "R",
    # Job was requeued.
    "REQUEUED": "RQ",
    # Job is about to change size.
    "RESIZING": "RS",
    # Sibling was removed from cluster due to other cluster starting the job.
    "REVOKED": "RV",
    # Job has an allocation, but execution has been suspended and CPUs have been
    # released for other jobs.
    "SUSPENDED": "S",
    # Job terminated upon reaching its time limit.
    "TIMEOUT": "TO",
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
        # whether the command failed. See the comment in call() if you want to
        # know more.
        call("scancel", "--verbose", job_id)

    def submit_target(self, target, dependencies):
        script = self.compile_script(target)
        args = ["--parsable"]
        if dependencies:
            args.append("--dependency=afterok:{}".format(":".join(dependencies)))

        environ = os.environ.copy()
        environ["GWF_EXEC_WORKFLOW_ROOT"] = self.working_dir
        return call("sbatch", *args, input=script, environ=environ).strip()

    def get_job_states_from_squeue(self, tracked_jobs):
        logger.debug("Loading job states from squeue")
        job_states = {}
        for line in call(
            "squeue", "--noheader", "--format=%i;%t", "--all"
        ).splitlines():
            job_id, state = line.strip().split(";")
            if job_id in tracked_jobs:
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
            # Slurm sometimes returns annoying stuff like "CANCELLED by 1234",
            # instead of a sensible, parsable value. So we do our best to clean
            # it up here.
            state = state.split()[0]
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
        job_states = {}
        if self.accounting_enabled:
            job_states.update(self.get_job_states_from_sacct_batched(tracked_jobs))
        job_states.update(self.get_job_states_from_squeue(tracked_jobs))
        return job_states

    def compile_script(self, target):
        with render_script(target) as buf:
            print(OPTION_STR.format("--job-name=", target.name), file=buf)
            for option_name, option_value in target.options.items():
                print(
                    OPTION_STR.format(OPTION_FLAGS[option_name], option_value), file=buf
                )

            stdout_flag = OPTION_STR.format(
                "--output=",
                os.path.join(self.working_dir, ".gwf", "logs", target.name + ".stdout"),
            )
            stderr_flag = OPTION_STR.format(
                "--error=",
                os.path.join(self.working_dir, ".gwf", "logs", target.name + ".stderr"),
            )
            devnull_flag = OPTION_STR.format("--output=", "/dev/null")
            if self.log_mode == "full":
                print(stdout_flag, file=buf)
                print(stderr_flag, file=buf)
            elif self.log_mode == "merged":
                print(stdout_flag, file=buf)
            elif self.log_mode == "none":
                print(devnull_flag, file=buf)

        return buf.getvalue()

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


@contextmanager
def render_script(target):
    buf = io.StringIO()
    print(f"#!{sys.executable} -mgwf.exec", file=buf)
    # print("# Generated by: gwf", file=buf)
    yield buf
    print(file=buf)
    serialize(target, buf)
