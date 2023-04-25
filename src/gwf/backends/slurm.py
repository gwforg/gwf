import logging
from collections import defaultdict

import attrs

from ..conf import config
from ..utils import ensure_trailing_newline, retry
from .base import PbsLikeBackendBase, BackendStatus
from .exceptions import BackendError
from .utils import call, has_exe

logger = logging.getLogger(__name__)


SLURM_RAW_STATES = {
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


SLURM_JOB_STATES = defaultdict(lambda: BackendStatus.UNKNOWN, SLURM_RAW_STATES)


@attrs.define
class SlurmBackend(PbsLikeBackendBase):
    """Backend for the Slurm workload manager.

    To use this backend you must activate the `slurm` backend.

    **Backend options:**

    * **backend.slurm.log_mode (str):** Must be either `full`, `merged` or
      `none`. If `full`, two log files will be stored for each target, one for
      standard output and one for standard error. If `merged`, only one log
      file will be written containing the combined streams. If `none`, no logs
      will be stored. (default: `full`).

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

    option_defaults = {
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

    option_flags = {
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

    option_str = "#SBATCH {0}{1}"

    @retry(on_exc=BackendError)
    def call_queue_command(self):
        return call("squeue", "--noheader", "--format=%i;%t", "--all")

    @retry(on_exc=BackendError)
    def call_cancel_command(self, job_id):
        # The --verbose flag here is necessary, otherwise we're not able to tell
        # whether the command failed. See the comment in call() if you
        # want to know more.
        return call("scancel", "--verbose", job_id)

    @retry(on_exc=BackendError)
    def call_submit_command(self, script, dependencies):
        args = ["--parsable"]
        if dependencies:
            args.append("--dependency=afterok:{}".format(":".join(dependencies)))
        return call("sbatch", *args, input=script)

    def parse_queue_output(self, stdout):
        job_states = {}
        for line in stdout.splitlines():
            job_id, state = line.split(";")
            job_states[job_id] = SLURM_JOB_STATES[state]
        return job_states

    def compile_script(self, target):
        out = []
        out.append("#!/bin/bash")
        out.append("# Generated by: gwf")

        out.append(self.option_str.format("--job-name=", target.name))

        for option_name, option_value in target.options.items():
            out.append(
                self.option_str.format(self.option_flags[option_name], option_value)
            )

        log_mode = config.get("backend.slurm.log_mode", "full")
        if log_mode == "full":
            out.append(
                self.option_str.format(
                    "--output=", self.log_manager.stdout_path(target)
                )
            )
            out.append(
                self.option_str.format("--error=", self.log_manager.stderr_path(target))
            )
        elif log_mode == "merged":
            out.append(
                self.option_str.format(
                    "--output=", self.log_manager.stdout_path(target)
                )
            )
        elif log_mode == "none":
            out.append(self.option_str.format("--output=", "/dev/null"))

        out.append("")
        out.append("cd {}".format(target.working_dir))
        out.append("export GWF_JOBID=$SLURM_JOBID")
        out.append('export GWF_TARGET_NAME="{}"'.format(target.name))
        out.append("set -e")
        out.append("")
        out.append(ensure_trailing_newline(target.spec))
        return "\n".join(out)

    @staticmethod
    def priority():
        if has_exe("sbatch") and has_exe("sinfo"):
            return 100
        return -100
