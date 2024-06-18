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

"""

import logging
import os.path
import re

import attrs

from ..utils import ensure_trailing_newline
from .base import BackendStatus, TrackingBackend
from .utils import call, has_exe

logger = logging.getLogger(__name__)

TARGET_DEFAULTS = {"queue": "normal", "memory": "4GB", "cores": 1}

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
        logger.debug(f"Cancelling job { job_id }")
        call("bkill", job_id)

    def submit_target(self, target, dependencies):
        script = self.compile_script(target)
        with open(
            os.path.join(self.working_dir, ".gwf", "logs", target.name + ".sh"), "w"
        ) as f:
            f.write(script)
        args = []
        if dependencies:
            args.append("-w")
            args.append(" && ".join([f"done({job_id})" for job_id in dependencies]))
        logger.debug(f"Submitting job { target.name } to LSF")
        stdout = call("bsub", *args, input=script).strip()
        job_id = re.search(r"Job <(\d+)>", stdout)[1]
        return job_id

    def get_job_states(self, tracked_jobs):
        logger.debug("Getting job states from LSF")
        if not tracked_jobs:
            return {}
        job_states = {job_id: BackendStatus.UNKNOWN for job_id in tracked_jobs}
        for job_id in tracked_jobs:
            cmd = ["bjobs", "-noheader", "-o", "stat", job_id]
            ret = call(*cmd).strip()
            if ret == "":
                continue
            state = BJOB_STATES[ret]
            job_states[job_id] = state
        return job_states

    def compile_script(self, target):
        target_options = target.options
        target_options["std_err"] = os.path.join(
            self.working_dir, ".gwf", "logs", target.name + ".stderr"
        )
        target_options["std_out"] = os.path.join(
            self.working_dir, ".gwf", "logs", target.name + ".stdout"
        )
        header = BJOB_HEADER
        for name, value in target_options.items():
            header = header.replace(f"{{{name}}}", str(value))
        header = header.replace("{job_name}", target.name)
        out = []
        out.append("#!/bin/bash")
        out.append(header)
        out.append("")
        out.append("# Generated by: gwf")
        out.append("")
        out.append("cd {}".format(target.working_dir))
        out.append("set -e")
        out.append("")
        out.append(ensure_trailing_newline(target.spec))

        return "\n".join(out)

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
