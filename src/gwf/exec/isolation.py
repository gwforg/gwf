import logging
import os
import shutil
import tempfile
import time
from collections.abc import Mapping
from contextlib import contextmanager
from pathlib import Path

import attrs

from gwf.exceptions import GWFError


logger = logging.getLogger(__name__)


class IsolationError(GWFError):
    pass


def _validate_mode(valid_modes):
    def validator(instance, attribute, value):
        if value not in valid_modes:
            modes = ", ".join(repr(mode) for mode in valid_modes)
            raise IsolationError(
                f"IsolationConfig {attribute.name} must be one of {modes}, "
                f"got {value!r}."
            )

    return validator


@attrs.frozen
class IsolationConfig:
    """Runtime configuration for an isolated target."""

    inputs: str = attrs.field(
        default="symlink", validator=_validate_mode(("copy", "symlink"))
    )
    outputs: str = attrs.field(
        default="move", validator=_validate_mode(("copy", "move"))
    )
    root: str | os.PathLike | None = attrs.field(
        default=None,
        validator=attrs.validators.optional(
            attrs.validators.instance_of((str, os.PathLike))
        ),
    )
    network: bool = attrs.field(default=False)


def isolate_network(command, target_name):
    """Run a command in a new user and network namespace."""
    logger.info("Target %s: blocking network access", target_name)
    return [
        "unshare",
        "--user",
        "--map-current-user",
        "--net",
        "--",
        *command,
    ]


def _flatten_paths(paths):
    if isinstance(paths, str) or hasattr(paths, "__fspath__"):
        yield paths
    elif isinstance(paths, Mapping):
        for value in paths.values():
            yield from _flatten_paths(value)
    else:
        for value in paths:
            yield from _flatten_paths(value)


def _relative_paths(paths, target_name):
    result = []
    seen = set()
    for declared_path in _flatten_paths(paths):
        path = Path(os.fspath(declared_path))
        if path.is_absolute():
            raise IsolationError(
                f"Target {target_name!r} cannot use absolute path {path!s} "
                "with isolation."
            )
        if not path.parts or path == Path(".") or ".." in path.parts:
            raise IsolationError(
                f"Target {target_name!r} cannot use path {path!s} with isolation."
            )

        path = Path(*(part for part in path.parts if part != "."))
        if path not in seen:
            result.append(path)
            seen.add(path)
    return result


def _validate_output_layout(paths, target_name):
    for index, path in enumerate(paths):
        for other in paths[index + 1 :]:
            if path in other.parents or other in path.parents:
                raise IsolationError(
                    f"Target {target_name!r} has overlapping isolated outputs "
                    f"{path!s} and {other!s}."
                )


def _validate_input_output_layout(inputs, outputs, target_name):
    for input_path in inputs:
        for output_path in outputs:
            if input_path in output_path.parents or output_path in input_path.parents:
                raise IsolationError(
                    f"Target {target_name!r} has overlapping isolated input "
                    f"{input_path!s} and output {output_path!s}."
                )


def _elapsed(start):
    return time.monotonic() - start


def _temporary_path(destination, marker):
    descriptor, path = tempfile.mkstemp(
        prefix=f".{destination.name}.gwf-{marker}-",
        dir=destination.parent,
    )
    os.close(descriptor)
    return Path(path)


def _unlink(path):
    try:
        path.unlink()
    except FileNotFoundError:
        pass


@attrs.define
class UnisolatedExecution:
    working_dir: Path
    isolated = False
    spec_directory = None

    def publish_outputs(self):
        pass


@attrs.define
class IsolatedExecution:
    target: object
    root: Path
    source_root: Path = attrs.field(init=False)
    working_dir: Path = attrs.field(init=False)
    spec_directory: Path = attrs.field(init=False)
    inputs: list = attrs.field(init=False)
    outputs: list = attrs.field(init=False)
    isolated = True

    def __attrs_post_init__(self):
        self.source_root = Path(self.target.working_dir).resolve()
        self.working_dir = self.root / "work"
        self.spec_directory = self.root
        self.inputs = _relative_paths(self.target.inputs, self.target.name)
        self.outputs = _relative_paths(self.target.outputs, self.target.name)
        _validate_output_layout(self.outputs, self.target.name)
        _validate_input_output_layout(self.inputs, self.outputs, self.target.name)

        if self.target.isolation.inputs == "symlink":
            shared_paths = set(self.inputs) & set(self.outputs)
            if shared_paths:
                paths = ", ".join(str(path) for path in sorted(shared_paths))
                raise IsolationError(
                    f"Target {self.target.name!r} cannot use symlinked inputs "
                    f"as outputs: {paths}."
                )

    def stage_inputs(self):
        self.working_dir.mkdir()
        logger.info(
            "Target %s: created isolated working directory %s",
            self.target.name,
            self.working_dir,
        )

        started = time.monotonic()
        for path in self.inputs:
            source = self.source_root / path
            destination = self.working_dir / path
            if not source.exists():
                raise IsolationError(
                    f"Input {path!s} for target {self.target.name!r} does not exist."
                )
            if not source.is_file():
                raise IsolationError(
                    f"Input {path!s} for target {self.target.name!r} is not a file."
                )

            destination.parent.mkdir(parents=True, exist_ok=True)
            transfer_started = time.monotonic()
            if self.target.isolation.inputs == "copy":
                logger.info(
                    "Target %s: copying input %s to %s",
                    self.target.name,
                    source,
                    destination,
                )
                shutil.copy2(source, destination)
                operation = "copied"
            else:
                logger.info(
                    "Target %s: symlinking input %s to %s",
                    self.target.name,
                    source,
                    destination,
                )
                destination.symlink_to(source)
                operation = "symlinked"
            logger.info(
                "Target %s: %s input %s to %s in %.3fs",
                self.target.name,
                operation,
                source,
                destination,
                _elapsed(transfer_started),
            )

        for path in self.outputs:
            (self.working_dir / path).parent.mkdir(parents=True, exist_ok=True)

        logger.info(
            "Target %s: staged %d input(s) in %.3fs",
            self.target.name,
            len(self.inputs),
            _elapsed(started),
        )

    def _validate_outputs(self):
        for path in self.outputs:
            source = self.working_dir / path
            if not source.exists():
                raise IsolationError(
                    f"Target {self.target.name!r} completed without producing "
                    f"declared output {path!s}."
                )
            if source.is_symlink() or not source.is_file():
                raise IsolationError(
                    f"Output {path!s} from target {self.target.name!r} is not "
                    "a regular file."
                )
            if not source.resolve().is_relative_to(self.working_dir.resolve()):
                raise IsolationError(
                    f"Output {path!s} from target {self.target.name!r} resolves "
                    "outside its isolated working directory."
                )

    def publish_outputs(self):
        self._validate_outputs()
        started = time.monotonic()
        for path in self.outputs:
            source = self.working_dir / path
            destination = self.source_root / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.is_dir() and not destination.is_symlink():
                raise IsolationError(
                    f"Cannot transfer output {path!s} from target "
                    f"{self.target.name!r} over a directory."
                )

            temporary_path = _temporary_path(destination, "output")
            transfer_started = time.monotonic()
            try:
                if self.target.isolation.outputs == "copy":
                    logger.info(
                        "Target %s: copying output %s to %s",
                        self.target.name,
                        source,
                        destination,
                    )
                    shutil.copy2(source, temporary_path)
                    operation = "copied"
                else:
                    logger.info(
                        "Target %s: moving output %s to %s",
                        self.target.name,
                        source,
                        destination,
                    )
                    shutil.move(source, temporary_path)
                    operation = "moved"
                os.replace(temporary_path, destination)
            finally:
                _unlink(temporary_path)

            logger.info(
                "Target %s: %s output %s to %s in %.3fs",
                self.target.name,
                operation,
                source,
                destination,
                _elapsed(transfer_started),
            )

        logger.info(
            "Target %s: transferred %d output(s) in %.3fs",
            self.target.name,
            len(self.outputs),
            _elapsed(started),
        )


@contextmanager
def execution_context(target):
    """Prepare and clean up the working directory for a target execution."""
    if target.isolation is None:
        yield UnisolatedExecution(Path(target.working_dir))
        return

    logger.setLevel(logging.INFO)
    temporary_root = target.isolation.root
    if temporary_root is not None:
        temporary_root = Path(temporary_root)
        if not temporary_root.is_absolute():
            temporary_root = Path(target.working_dir) / temporary_root

    try:
        temporary_directory = tempfile.TemporaryDirectory(
            prefix=f"gwf-{target.name}-",
            dir=temporary_root,
        )
    except OSError as error:
        location = temporary_root or "the system temporary directory"
        raise IsolationError(
            f"Could not create an isolated working directory for target "
            f"{target.name!r} under {location}: {error}"
        ) from error

    with temporary_directory as root:
        execution = IsolatedExecution(target, Path(root))
        execution.stage_inputs()
        try:
            yield execution
        finally:
            logger.info(
                "Target %s: removing isolated working directory %s",
                target.name,
                execution.working_dir,
            )
