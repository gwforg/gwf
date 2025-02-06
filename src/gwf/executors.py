import attrs
import base64
import logging
import os
import pickle
import prettyprinter
import shutil
from typing import Optional, Protocol, Iterable

from .exceptions import GWFError


prettyprinter.install_extras(["attrs"])

logger = logging.getLogger(__name__)


class Executor(Protocol):
    def get_command(self, spec_path: str, workflow_root: str) -> Iterable[str]:
        return []


@attrs.define
class Bash:
    def get_command(self, spec_path: str, workflow_root: str) -> Iterable[str]:
        return [spec_path]


@attrs.define
class Conda:
    env: str = attrs.field()
    debug_mode: bool = attrs.field(default=False)

    def get_command(self, spec_path: str, workflow_root: str) -> Iterable[str]:
        debug_flags = ["--debug-wrapper-scripts", "-vvv"] if self.debug_mode else []
        conda_exe = os.environ.get("CONDA_EXE", shutil.which("conda"))
        if conda_exe is None:
            raise GWFError("Could not find conda installation")
        logger.debug("found conda executable at %s", conda_exe)
        return [conda_exe, "run", "--live-stream", *debug_flags, "-n", self.env, spec_path]


@attrs.define
class Pixi:
    project: Optional[str] = attrs.field(default=None)
    env: str = attrs.field(default="default")
    debug_mode: bool = attrs.field(default=False)

    def get_command(self, spec_path: str, workflow_root: str) -> Iterable[str]:
        pixi_exe = os.getenv("PIXI_EXE", shutil.which("pixi"))
        if pixi_exe is None:
            raise GWFError("Could not find pixi installation")
        debug_flags = ["-vvv"] if self.debug_mode else []
        manifest_path = workflow_root if self.project is None else self.project
        return [pixi_exe, "run", "--no-install", "--frozen", *debug_flags, "-e", self.env, "--manifest-path", manifest_path, spec_path]


@attrs.define
class Singularity:
    image: str = attrs.field()
    flags: Iterable[str] = attrs.field(factory=list)
    debug_mode: bool = attrs.field(default=False)

    def get_command(self, spec_path: str, workflow_root: str) -> Iterable[str]:
        singularity_exe = shutil.which("singularity")
        if singularity_exe is None:
            raise GWFError("Could not find singularity installation")
        debug_flags = ["--debug"] if self.debug_mode else []
        return [singularity_exe, *debug_flags, "exec", *self.flags, self.image, spec_path]


@attrs.define
class Apptainer:
    image: str = attrs.field()
    flags: Iterable[str] = attrs.field(factory=list)
    debug_mode: bool = attrs.field(default=False)

    def get_command(self, spec_path: str, workflow_root: str) -> Iterable[str]:
        apptainer_exe = shutil.which("apptainer")
        if apptainer_exe is None:
            raise GWFError("Could not find apptainer installation")
        debug_flags = ["--debug"] if self.debug_mode else []
        return [apptainer_exe, *debug_flags, "exec", *self.flags, self.image, spec_path]


def skip_comments(script_file):
    for line in script_file:
        if line.startswith("#GWF COMMENT"):
            continue
        yield line


def consume(cond, lst):
    i = 0
    while cond(lst[i]):
        i += 1
    return lst[i:]


def take(cond, lst):
    i = 0
    buf = []
    while cond(lst[i]):
        buf.append(lst[i])
        i += 1
    return buf, lst[i:]


def serialize(target, file):
    print("#GWF TARGET", file=file)
    repr_str = prettyprinter.pformat(target)
    for line in repr_str.splitlines(keepends=False):
        print(f"#GWF COMMENT {line}", file=file)
    print("#GWF SPEC", file=file)
    print(file=file)
    print(target.spec, file=file)
    print(file=file)

    data_width = 70
    pickled = base64.b64encode(pickle.dumps(target)).decode("ascii")
    for i in range(0, len(pickled), data_width):
        print("#GWF DATA", pickled[i:i + data_width], file=file)
    print("#GWF END", file=file)
    print(file=file)


def deserialize(script_file):
    script_file = list(skip_comments(script_file))
    script_file = consume(lambda line: not line.startswith("#GWF TARGET"), script_file)
    script_file = consume(lambda line: line.startswith("#GWF TARGET"), script_file)
    script_file = consume(lambda line: line.startswith("#GWF SPEC"), script_file)
    script_file = consume(lambda line: not line.startswith("#GWF DATA"), script_file)

    data_lines, script_file = take(lambda line: line.startswith("#GWF DATA"), script_file)
    data = "".join(line.replace("#GWF DATA ", "") for line in data_lines).replace("\n", "")

    script_file = consume(lambda line: line.startswith("#GWF END"), script_file)
    return pickle.loads(base64.b64decode(data))
