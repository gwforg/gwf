"""Backend for Sun Grid Engine (SGE).

.. versionchanged:: 3.0.0
    ``sge_args`` was added for passing extra arguments directly to
    :command:`qsub`, and executors are now supported.

To use this backend you must activate the `sge` backend. The backend
currently assumes that a SGE parallel environment called "smp" is
available. You can check which parallel environments are available on your
system by running :command:`qconf -spl`.

.. note::
    Dependencies are expressed with :command:`qsub -hold_jid`, which only
    waits for the dependencies to *finish*, not to succeed. Unlike the other
    backends, a target may therefore start even though one of its
    dependencies failed, and will then usually fail itself because of missing
    input files. When several targets have failed, look for the failed target
    furthest upstream.

**Backend options:**

None.

**Target options:**

* **cores (int):**
    Number of cores allocated to this target (default: 1).
* **memory (str):**
    Total memory allocated to this target (default: 1g). *gwf* divides it
    between the target's cores, since SGE expects memory per core.
* **walltime (str):**
    Time limit for this target (default: 01:00:00).
* **queue (str):**
    Queue to submit the target to. To specify multiple queues, specify a
    comma-separated list of queue names.
* **account (str):**
    Account to be used when running the target. Corresponds to the SGE
    project.
* **sge_args (list[str]):**
    Extra arguments passed directly to `qsub`. Use this for SGE options that
    are not otherwise supported by *gwf*.
"""

import logging
import math
import os
import re
from decimal import Decimal
from xml.etree import ElementTree

import attrs

from ..executors import render_exec_script
from ..log_storage import get_log_paths
from .base import BackendStatus, TrackingBackend
from .exceptions import BackendError
from .utils import call, has_exe

logger = logging.getLogger(__name__)


TARGET_DEFAULTS = {
    "cores": 1,
    "memory": "1g",
    "walltime": "01:00:00",
    "queue": None,
    "account": None,
    "sge_args": [],
}

SUBMISSION_ARGS_OPTION = "sge_args"

OPTION_FLAGS = {
    "cores": "-pe smp ",
    "memory": "-l h_vmem=",
    "walltime": "-l h_rt=",
    "queue": "-q ",
    "account": "-P ",
}


MEMORY_PATTERN = re.compile(r"(\d+(?:\.\d+)?)\s*([kmgKMG]?)[bB]?")

# In SGE, lowercase memory units are powers of 1000 and uppercase units are
# powers of 1024.
MEMORY_UNITS = ("", "k", "m", "g")


def per_core_memory(memory, cores):
    """Return total `memory` divided between `cores` as an SGE memory value.

    The result keeps the unit of `memory` when the division is exact.
    Otherwise it moves to smaller units, rounding up to whole bytes so that the
    job is never given less than the total memory requested.
    """
    match = MEMORY_PATTERN.fullmatch(str(memory).strip())
    if match is None:
        raise BackendError(f"Invalid memory value: {memory!r}")
    amount, unit = match.groups()
    base = 1024 if unit.isupper() else 1000

    per_core = Decimal(amount) / cores
    index = MEMORY_UNITS.index(unit.lower())
    while index > 0 and per_core != per_core.to_integral_value():
        per_core *= base
        index -= 1

    unit = MEMORY_UNITS[index]
    if base == 1024:
        unit = unit.upper()
    return f"{math.ceil(per_core)}{unit}"


@attrs.define
class SGEOps:
    working_dir: str = attrs.field()
    target_defaults: dict = attrs.field()

    def cancel_job(self, job_id):
        # The --verbose flag here is necessary, otherwise we're not able to tell
        # whether the command failed. See the comment in call() if you
        # want to know more.
        return call("qdel", job_id)

    def submit_target(self, target, dependencies):
        script = self.compile_script(target)
        args = ["-terse"]
        if dependencies:
            args.append("-hold_jid")
            args.append(",".join(dependencies))
        args.extend(target.options.get(SUBMISSION_ARGS_OPTION, []))
        environ = os.environ.copy()
        environ["GWF_EXEC_WORKFLOW_ROOT"] = self.working_dir
        # `qsub -terse` prints the job id followed by a newline.
        return call("qsub", *args, input=script, environ=environ).strip()

    def get_job_states(self, tracked_jobs):
        job_states = {}
        root = ElementTree.fromstring(call("qstat", "-f", "-xml"))
        for job in root.iter("job_list"):
            job_id = job.find("JB_job_number").text
            state = job.find("state").text
            assert job_id is not None
            assert state is not None

            # Guessing job state based on
            # https://gist.github.com/cmaureir/4fa2d34bc9a1bd194af1
            if "d" in state or "E" in state:
                job_state = BackendStatus.UNKNOWN
            elif "r" in state or "t" in state or "s" in state:
                job_state = BackendStatus.RUNNING
            else:
                job_state = BackendStatus.SUBMITTED
            job_states[job_id] = job_state
        return job_states

    def compile_script(self, target):
        option_str = "#$ {0}{1}"
        with render_exec_script(target) as buf:
            print(option_str.format("-N ", target.name), file=buf)
            print("#$ -V", file=buf)
            print("#$ -w v", file=buf)
            print("#$ -cwd", file=buf)

            for option_name, option_value in target.options.items():
                if option_name == SUBMISSION_ARGS_OPTION:
                    continue
                # SGE wants per-core memory, but gwf wants total memory.
                if option_name == "memory":
                    option_value = per_core_memory(
                        option_value, target.options.get("cores", 1)
                    )
                print(option_str.format(OPTION_FLAGS[option_name], option_value), file=buf)

            stdout_path, stderr_path = get_log_paths(self.working_dir, target.name)
            print(option_str.format("-o ", stdout_path), file=buf)
            print(option_str.format("-e ", stderr_path), file=buf)

        return buf.getvalue()

    def close(self):
        pass


def create_backend(working_dir):
    return TrackingBackend(
        working_dir,
        name="sge",
        ops=SGEOps(working_dir, target_defaults=TARGET_DEFAULTS),
    )


def priority():
    if has_exe("qsub") and has_exe("qdel"):
        return 50
    return -100


setup = (create_backend, priority())
