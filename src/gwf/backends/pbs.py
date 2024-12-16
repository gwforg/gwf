"""
Backend for Portable Batch System (PBS).

To use this backend, you must activate the `pbs` backend. This backend requires
the commands `qsub` and `qstat`.

**Backend options:**

None.

**Target options:**

* **cores (int):**
    Number of cores allocated to this target (default: 1).
* **memory (str):**
    Memory allocated to this target (default: 4GB).
* **queue (str):**
    Queue to submit the target to (default: normal).

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

PBS_HEADER = """#PBS -l mem={memory}
#PBS -l nodes=1:ppn={cores}
#PBS -q {queue}
#PBS -o {std_out}
#PBS -e {std_err}
#PBS -N {job_name}"""

PBS_STATES = {
    "Q": BackendStatus.SUBMITTED,
    "H": BackendStatus.SUBMITTED,
    "R": BackendStatus.RUNNING,
    "E": BackendStatus.RUNNING,
    "C": BackendStatus.COMPLETED,
    "F": BackendStatus.FAILED,
    "S": BackendStatus.FAILED,
    "U": BackendStatus.UNKNOWN,
}


@attrs.define
class PBSOps:
    working_dir: str = attrs.field()
    target_defaults: dict = attrs.field()

    def cancel_job(self, job_id):
        logger.debug(f"Cancelling job { job_id }")
        call("qdel", job_id)

    def submit_target(self, target, dependencies):
        script = self.compile_script(target)
        script_path = os.path.join(
            self.working_dir, ".gwf", "logs", target.name + ".sh"
        )
        with open(script_path, "w") as f:
            f.write(script)
        args = []
        if dependencies:
            args.append("-W depend=afterok:" + ":".join(dependencies))
        logger.debug(f"Submitting job { target.name } to PBS")
        stdout = call("qsub", *args, script_path).strip()
        job_id = stdout.split(".")[0]  # Extract job ID from the full PBS job ID
        return job_id

    def get_job_states(self, tracked_jobs):
        logger.debug("Getting job states from PBS")
        if not tracked_jobs:
            return {}
        job_states = {job_id: BackendStatus.UNKNOWN for job_id in tracked_jobs}
        for job_id in tracked_jobs:
            cmd = ["qstat", "-f", job_id]
            ret = call(*cmd)
            if "job_state" not in ret:
                continue
            state = re.search(r"job_state = (\w)", ret).group(1)
            job_states[job_id] = PBS_STATES.get(state, BackendStatus.UNKNOWN)
        return job_states

    def compile_script(self, target):
        target_options = target.options
        target_options["std_err"] = os.path.join(
            self.working_dir, ".gwf", "logs", target.name + ".stderr"
        )
        target_options["std_out"] = os.path.join(
            self.working_dir, ".gwf", "logs", target.name + ".stdout"
        )
        header = PBS_HEADER
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
        name="pbs",
        ops=PBSOps(working_dir, target_defaults=TARGET_DEFAULTS),
    )


def priority():
    if has_exe("qsub"):
        return 100
    return -100


setup = (create_backend, priority())