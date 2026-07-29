"""Backend for IBM Load Sharing Facility (LSF).

To use this backend, you must activate the `lsf` backend. This backend requires
the commands `bsub` and `bjobs`.

**Backend options:**

None.

**Target options:**

* **cores (int):**
    Number of cores allocated to this target (default: 1).
* **memory (str):**
    Memory allocated to this target (default: 4GB).
* **queue (str):**
    Queue to submit the target to (default: normal).
    for different purposes or priorities.
* **lsf_args (list[str]):**
    Extra arguments passed directly to `bsub`. Use this for LSF options that
    are not otherwise supported by *gwf*.

"""

import logging
import os
import re

import attrs

from ..executors import render_exec_script
from ..log_storage import get_log_paths
from .base import BackendStatus, TrackingBackend
from .utils import call, has_exe

logger = logging.getLogger(__name__)

TARGET_DEFAULTS = {
    "queue": "normal",
    "memory": "4GB",
    "cores": 1,
    "lsf_args": [],
}

SUBMISSION_ARGS_OPTION = "lsf_args"

BJOB_HEADER = """#BSUB -M {memory}
#BSUB -R "select[mem>{memory}] rusage[mem={memory}] span[hosts=1]"
#BSUB -n {cores}
#BSUB -q {queue}
#BSUB -oo {std_out}
#BSUB -eo {std_err}
#BSUB -J {job_name}"""

BJOB_STATES = {
    "PEND": BackendStatus.SUBMITTED,
    "WAIT": BackendStatus.SUBMITTED,
    "RUN": BackendStatus.RUNNING,
    "ZOMBI": BackendStatus.RUNNING,
    "DONE": BackendStatus.COMPLETED,
    "EXIT": BackendStatus.FAILED,
    "PSUSP": BackendStatus.FAILED,
    "USUSP": BackendStatus.FAILED,
    "SSUSP": BackendStatus.FAILED,
    "UNKWN": BackendStatus.UNKNOWN,
}


@attrs.define
class LSFOps:
    working_dir: str = attrs.field()
    target_defaults: dict = attrs.field()

    def cancel_job(self, job_id):
        logger.debug(f"Cancelling job {job_id}")
        call("bkill", job_id)

    def submit_target(self, target, dependencies):
        script = self.compile_script(target)
        args = []
        if dependencies:
            args.append("-w")
            args.append(" && ".join([f"done({job_id})" for job_id in dependencies]))
            args.append("-ti")
        args.extend(target.options.get(SUBMISSION_ARGS_OPTION, []))
        logger.debug(f"Submitting job {target.name} to LSF")
        environ = os.environ.copy()
        environ["GWF_EXEC_WORKFLOW_ROOT"] = self.working_dir
        stdout = call("bsub", *args, input=script, environ=environ).strip()
        job_id = re.search(r"Job <(\d+)>", stdout)[1]
        return job_id

    def get_job_states(self, tracked_jobs):
        logger.debug("Getting job states from LSF")
        if not tracked_jobs:
            return {}
        # Default tracked job IDs to UNKNOWN
        job_states = {job_id: BackendStatus.UNKNOWN for job_id in tracked_jobs}

        # Get all current job statuses
        for line in call("bjobs", "-noheader", "-o" , "jobid stat delimiter=','").splitlines():
            job_id, state = line.strip().split(",")
            # Update job status if tracked
            if job_id in tracked_jobs:
                job_states[job_id] = BJOB_STATES[state]
        return job_states

    def compile_script(self, target):
        target_options = target.options
        stdout_path, stderr_path = get_log_paths(self.working_dir, target.name)
        target_options["std_err"] = stderr_path
        target_options["std_out"] = stdout_path
        header = BJOB_HEADER
        for name, value in target_options.items():
            header = header.replace(f"{{{name}}}", str(value))
        header = header.replace("{job_name}", target.name)
        with render_exec_script(target) as buf:
            print(header, file=buf)
        return buf.getvalue()

    def close(self):
        pass


def create_backend(working_dir):
    return TrackingBackend(
        working_dir,
        name="lsf",
        ops=LSFOps(working_dir, target_defaults=TARGET_DEFAULTS),
    )


def priority():
    if has_exe("bsub"):
        return 100
    return -100


setup = (create_backend, priority())
