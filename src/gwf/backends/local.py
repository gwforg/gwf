"""Backend that runs targets on a local cluster.

To use this backend you must activate the `local` backend and start a local
scheduler that the backend can submit targets to. To make two cores available
to the scheduler, run the command::

    gwf -b local workers -n 2

in the working directory of your project. The workflow file must be accessible
to *gwf*. Thus, if your workflow file is not called `workflow.py` or the
workflow object is not called `gwf`, you must specify this so that *gwf* can
locate the workflow::

    gwf -f myworkflow.py:wf1 -b local workers -n 2

If the local backend is your default backend you can of course omit the
``-b local`` option.

If the ``-n`` option is omitted, *gwf* will detect the number of cores available
and use all of them. The scheduler uses the system's physical memory by default;
use ``--max-memory`` to set another capacity.

To run your workflow, open another terminal and then type::

    gwf -b local run

To stop the scheduler press :kbd:`Control-c`.

The scheduler communicates over a per-workflow Unix domain socket under
``$XDG_RUNTIME_DIR/gwf``. If ``XDG_RUNTIME_DIR`` is unavailable, *gwf* uses a
private directory under the system temporary directory.

**Target options:**

* **cores (int):**
    Number of cores allocated to this target (default: 1).
* **memory (str):**
    Memory allocated to this target (default: 1g).
* **walltime (str):**
    Wall-clock time limit, ``HH:MM:SS`` (default: ``01:00:00``).
"""

import asyncio
from contextlib import contextmanager, suppress
import fcntl
import hashlib
import itertools
import io
import json
import logging
import multiprocessing
import os
import re
import signal
import socket
import sys
import tempfile
from enum import Enum
from io import TextIOWrapper
from pathlib import Path

import attrs

from ..log_storage import get_log_paths
from ..executors import serialize
from .base import BackendStatus, TrackingBackend
from .exceptions import BackendError

__all__ = ("Client", "Server")


logger = logging.getLogger(__name__)

TARGET_DEFAULTS = {
    "cores": 1,
    "memory": "1g",
    "walltime": "01:00:00",
}

MEMORY_UNITS = {
    "": 1,
    "k": 1024,
    "m": 1024**2,
    "g": 1024**3,
    "t": 1024**4,
}

MEMORY_PATTERN = re.compile(r"(\d+(?:\.\d+)?)\s*([kmgt]?)(?:i?b)?", re.IGNORECASE)


def parse_memory(value):
    match = MEMORY_PATTERN.fullmatch(str(value).strip())
    if match is None:
        raise ValueError(f"Invalid memory value: {value!r}")
    amount, unit = match.groups()
    return int(float(amount) * MEMORY_UNITS[unit.lower()])


def format_memory(value):
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.1f}{unit}"
        value /= 1024


def get_total_memory():
    return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")


def parse_max_memory(value):
    if value is None:
        return get_total_memory()
    return parse_memory(value)


def parse_walltime(value):
    if value is None:
        return None
    fields = str(value).split(":")
    if len(fields) == 1:
        return float(fields[0])
    hours, minutes, seconds = fields
    return int(hours) * 60 * 60 + int(minutes) * 60 + float(seconds)


def _get_runtime_directory():
    runtime_directory = os.environ.get("XDG_RUNTIME_DIR")
    if runtime_directory:
        app_runtime_directory = Path(runtime_directory) / "gwf"
    else:
        app_runtime_directory = (
            Path(tempfile.gettempdir()).resolve() / f"gwf-{os.geteuid()}"
        )
    app_runtime_directory.mkdir(mode=0o700, parents=True, exist_ok=True)
    return app_runtime_directory


def get_socket_path(working_dir):
    workflow_root = Path(working_dir).resolve()
    workflow_hash = hashlib.sha256(os.fsencode(workflow_root)).hexdigest()[:16]
    return _get_runtime_directory() / f"{workflow_hash}.sock"


@contextmanager
def _socket_lock(socket_path):
    lock_path = socket_path.with_suffix(".lock")
    with lock_path.open("w") as lock_file:
        try:
            fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise BackendError(
                f"A local scheduler is already running at {socket_path}."
            ) from error
        yield


class LocalStatus(Enum):
    UNKNOWN = 0  # task is unknown
    SUBMITTED = 1  # task was submitted
    RUNNING = 2  # task is running
    FAILED = 3  # script failed with non-zero return code
    COMPLETED = 4  # completed successfully
    CANCELLED = 5  # cancelled by user
    KILLED = 6  # killed because of timeout


STATUS_MAP = {
    LocalStatus.UNKNOWN: BackendStatus.UNKNOWN,
    LocalStatus.SUBMITTED: BackendStatus.SUBMITTED,
    LocalStatus.RUNNING: BackendStatus.RUNNING,
    LocalStatus.FAILED: BackendStatus.FAILED,
    LocalStatus.COMPLETED: BackendStatus.COMPLETED,
    LocalStatus.CANCELLED: BackendStatus.CANCELLED,
    LocalStatus.KILLED: BackendStatus.FAILED,
}


class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, LocalStatus):
            return obj.name
        return json.JSONEncoder.default(self, obj)


def decode(data):
    msg = json.loads(data.strip())
    kind = msg.pop("__kind__")
    return kind, msg


def encode(kind, **kwargs):
    payload = dict(__kind__=kind, **kwargs)
    data = json.dumps(payload, cls=CustomEncoder) + "\n"
    return data


@attrs.frozen
class Client:
    """A client for communicating with the local backend server."""

    sock: socket.socket = attrs.field(hash=False)
    reader: TextIOWrapper = attrs.field(repr=False)
    writer: TextIOWrapper = attrs.field(repr=False)

    @classmethod
    def from_socket(cls, sock):
        reader = sock.makefile(encoding="utf-8", mode="r")
        writer = sock.makefile(encoding="utf-8", mode="w")
        return cls(sock=sock, reader=reader, writer=writer)

    def send(self, kind, **msg):
        self.writer.write(encode(kind, **msg))
        self.writer.flush()

    def recv(self):
        return decode(self.reader.readline())

    @classmethod
    def connect(cls, socket_path):
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        try:
            sock.connect(os.fspath(socket_path))
        except OSError:
            sock.close()
            raise
        return cls.from_socket(sock)

    def submit(self, target, deps=None):
        script = io.StringIO()
        serialize(target, script)
        self.send(
            "enqueue_task",
            name=target.name,
            script=script.getvalue(),
            time_limit=target.options.get("walltime", TARGET_DEFAULTS["walltime"]),
            working_dir=target.working_dir,
            deps=deps or [],
            cores=target.options.get("cores", TARGET_DEFAULTS["cores"]),
            memory=target.options.get("memory", TARGET_DEFAULTS["memory"]),
        )
        kind, response = self.recv()
        if kind == "task_rejected":
            raise BackendError(response["message"])
        assert kind == "task_enqueued"
        return response["tid"]

    def status(self):
        self.send("get_task_states")
        msg_type, response = self.recv()
        assert msg_type == "task_states", "invalid response received"
        return {k: LocalStatus[v] for k, v in response["tasks"].items()}

    def cancel(self, job_id):
        self.send("cancel_task", tid=job_id)
        msg_type, response = self.recv()
        assert msg_type == "task_cancelled", "invalid response received"
        assert response["tid"] == job_id

    def shutdown(self):
        self.send("shutdown")
        msg_type, response = self.recv()
        assert msg_type == "shutdown", "invalid response received"
        assert not response

    def close(self):
        with suppress(OSError):
            self.send("close")
        with suppress(OSError):
            self.writer.close()
        with suppress(OSError):
            self.reader.close()
        self.sock.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


@attrs.define
class LocalOps:
    working_dir: str = attrs.field()
    socket_path: Path = attrs.field(converter=Path)
    target_defaults: dict = attrs.field()

    _client: Client = attrs.field(init=False, repr=False)

    @_client.default
    def _create_client(self):
        try:
            return Client.connect(self.socket_path)
        except (ConnectionRefusedError, FileNotFoundError) as error:
            raise BackendError(
                "Local backend could not connect to the scheduler for "
                f"{self.working_dir} at {self.socket_path}."
            ) from error

    def get_job_states(self, tracked_jobs):
        return {
            int(k): STATUS_MAP[v]
            for k, v in self._client.status().items()
            if int(k) in tracked_jobs
        }

    def submit_target(self, target, dependency_ids):
        return self._client.submit(target, deps=dependency_ids)

    def cancel_job(self, job_id):
        self._client.cancel(job_id)

    def close(self):
        self._client.close()


def create_backend(working_dir):
    return TrackingBackend(
        working_dir,
        name="local",
        ops=LocalOps(
            working_dir,
            get_socket_path(working_dir),
            target_defaults=TARGET_DEFAULTS,
        ),
    )


@attrs.frozen
class TaskFailedError(Exception):
    returncode: int


@attrs.frozen
class TimeLimitExceededError(Exception):
    pass


class TargetLogger(logging.LoggerAdapter):
    def process(self, message, kwargs):
        prefix = f"Target {self.extra['name']} (id {self.extra['tid']})"
        return f"{prefix}: {message}", kwargs


@attrs.define
class ScheduledTarget:
    tid: int = attrs.field()
    name: str = attrs.field()
    script: str = attrs.field()
    working_dir: str = attrs.field()
    time_limit: float | None = attrs.field()
    deps: tuple = attrs.field(converter=tuple)
    cores: int = attrs.field()
    memory: int = attrs.field()
    done: asyncio.Future = attrs.field(repr=False)
    execution: asyncio.Task | None = attrs.field(default=None, repr=False)
    cancel_requested: bool = attrs.field(default=False)
    log: TargetLogger = attrs.field(init=False, repr=False)

    @log.default
    def _create_logger(self):
        return TargetLogger(logger, {"name": self.name, "tid": self.tid})


@attrs.define
class Scheduler:
    working_dir: Path = attrs.field(converter=Path)
    max_cores: int = attrs.field(default=multiprocessing.cpu_count())
    max_memory: int = attrs.field(default=None, converter=parse_max_memory)

    tid_generator: object = attrs.field(factory=itertools.count)
    task_states: dict = attrs.field(factory=dict)
    tasks: dict = attrs.field(factory=dict)
    accepting_tasks: bool = attrs.field(default=True)
    available_cores: int = attrs.field(init=False)
    available_memory: int = attrs.field(init=False)

    @available_cores.default
    def _init_available_cores(self):
        return self.max_cores

    @available_memory.default
    def _init_available_memory(self):
        return self.max_memory

    def enqueue_task(
        self,
        name,
        script,
        working_dir,
        time_limit,
        deps,
        cores=1,
        memory="1g",
    ):
        if not self.accepting_tasks:
            raise BackendError("Local scheduler is shutting down.")

        memory = parse_memory(memory)
        time_limit = parse_walltime(time_limit)

        if cores < 1:
            raise BackendError(
                f"Target {name} requested an invalid number of cores: {cores!r}."
            )
        if cores > self.max_cores:
            raise BackendError(
                f"Target {name} requested {cores} cores, but the local "
                f"scheduler has only {self.max_cores}."
            )
        if memory < 1:
            raise BackendError(
                f"Target {name} requested an invalid amount of memory: {memory}."
            )
        if memory > self.max_memory:
            raise BackendError(
                f"Target {name} requested {format_memory(memory)} of memory, "
                f"but the local scheduler has only {format_memory(self.max_memory)}."
            )

        tid = next(self.tid_generator)
        target = ScheduledTarget(
            tid=tid,
            name=name,
            script=script,
            working_dir=working_dir,
            time_limit=time_limit,
            deps=deps,
            cores=cores,
            memory=memory,
            done=asyncio.get_running_loop().create_future(),
        )
        self.tasks[tid] = target
        self.task_states[tid] = LocalStatus.SUBMITTED
        target.log.debug(
            "queued with %d core(s), %s memory, time limit %s, "
            "dependencies %s, working directory %s",
            target.cores,
            format_memory(target.memory),
            target.time_limit,
            list(target.deps),
            target.working_dir,
        )
        self._schedule()
        return tid

    async def cancel_task(self, tid):
        target = self.tasks[tid]
        state = self.task_states[tid]
        if state == LocalStatus.SUBMITTED:
            target.log.debug("cancelling before execution")
            target.cancel_requested = True
            self._finish_target(target, LocalStatus.CANCELLED)
            self._schedule()
            return
        if state != LocalStatus.RUNNING:
            target.log.debug(
                "ignoring cancellation request in state %s",
                state.name.lower(),
            )
            return

        if not target.cancel_requested:
            target.log.debug("cancelling running execution")
            target.cancel_requested = True
            target.execution.cancel()
        try:
            await asyncio.shield(target.done)
        except asyncio.CancelledError:
            await target.done
            raise

    async def kill(self):
        await asyncio.gather(
            *(self.cancel_task(tid) for tid in self.tasks),
            return_exceptions=True,
        )

    async def wait(self):
        if self.tasks:
            await asyncio.wait({target.done for target in self.tasks.values()})

    async def shutdown(self):
        self.accepting_tasks = False
        await self.kill()
        await self.wait()

    def get_task_state(self, tid):
        return self.task_states.get(tid, None)

    def get_task_states(self):
        return dict(self.task_states)

    async def wait_for(self, tids, timeout=None):
        done = {self.tasks[tid].done for tid in tids}
        return await asyncio.wait(done, timeout=timeout)

    def _finish_target(self, target, state):
        if target.done.done():
            return
        previous_state = self.task_states[target.tid]
        self.task_states[target.tid] = state
        target.done.set_result(state)
        target.log.debug(
            "state changed from %s to %s",
            previous_state.name.lower(),
            state.name.lower(),
        )

    def _schedule(self):
        for target in self.tasks.values():
            if self.task_states[target.tid] != LocalStatus.SUBMITTED:
                continue

            try:
                dependency_states = [
                    self.task_states[dep_tid] for dep_tid in target.deps
                ]
            except KeyError:
                target.log.exception(
                    "cannot schedule because dependencies %s include an unknown target",
                    list(target.deps),
                )
                self._finish_target(target, LocalStatus.FAILED)
                continue

            if any(
                state in (LocalStatus.SUBMITTED, LocalStatus.RUNNING)
                for state in dependency_states
            ):
                continue

            failed_dependency = next(
                (
                    (dep_tid, state)
                    for dep_tid, state in zip(target.deps, dependency_states)
                    if state != LocalStatus.COMPLETED
                ),
                None,
            )
            if failed_dependency is not None:
                dep_tid, dep_state = failed_dependency
                target.log.debug(
                    "will not run because dependency %d finished in state %s",
                    dep_tid,
                    dep_state.name.lower(),
                )
                self._finish_target(target, dep_state)
                continue

            if (
                target.cores > self.available_cores
                or target.memory > self.available_memory
            ):
                continue

            self.available_cores -= target.cores
            self.available_memory -= target.memory
            self.task_states[target.tid] = LocalStatus.RUNNING
            target.log.debug(
                "starting with %d core(s) and %s memory; "
                "%d of %d core(s) and %s of %s memory remain available",
                target.cores,
                format_memory(target.memory),
                self.available_cores,
                self.max_cores,
                format_memory(self.available_memory),
                format_memory(self.max_memory),
            )
            target.execution = asyncio.create_task(self._run_target(target))

    async def _run_target(self, target):
        try:
            await self._execute_task(target)
        except asyncio.CancelledError:
            target.log.debug("execution was cancelled")
            state = LocalStatus.CANCELLED
        except TimeLimitExceededError:
            target.log.debug("execution exceeded its time limit %s", target.time_limit)
            state = LocalStatus.KILLED
        except TaskFailedError as error:
            target.log.debug(
                "execution failed with return code %d",
                error.returncode,
            )
            state = LocalStatus.FAILED
        except Exception:
            target.log.exception("execution failed unexpectedly")
            state = LocalStatus.FAILED
        else:
            state = LocalStatus.COMPLETED
        finally:
            target.execution = None
            self.available_cores += target.cores
            self.available_memory += target.memory
            target.log.debug(
                "released %d core(s) and %s memory; "
                "%d of %d core(s) and %s of %s memory are available",
                target.cores,
                format_memory(target.memory),
                self.available_cores,
                self.max_cores,
                format_memory(self.available_memory),
                format_memory(self.max_memory),
            )

        self._finish_target(target, state)
        self._schedule()

    async def _terminate_process_group(
        self,
        proc,
        communication_task,
        target_log,
        grace_period=10,
    ):
        process_group_id = proc.pid
        target_log.debug("sending SIGTERM to process group %d", process_group_id)
        try:
            os.killpg(process_group_id, signal.SIGTERM)
        except ProcessLookupError:
            target_log.debug("process group %d has already exited", process_group_id)

        try:
            await asyncio.wait_for(
                asyncio.shield(communication_task),
                timeout=grace_period,
            )
        except asyncio.TimeoutError:
            target_log.debug(
                "process group %d did not exit after %d seconds",
                process_group_id,
                grace_period,
            )

        target_log.debug("sending SIGKILL to process group %d", process_group_id)
        try:
            os.killpg(process_group_id, signal.SIGKILL)
        except PermissionError:
            target_log.debug(
                "permission denied when signalling process group %d",
                process_group_id,
            )
        except ProcessLookupError:
            target_log.debug("process group %d has exited", process_group_id)

        try:
            await communication_task
        except Exception:
            target_log.debug("failed while draining process output", exc_info=True)
        await proc.wait()
        target_log.debug(
            "process %d exited with return code %s",
            proc.pid,
            proc.returncode,
        )

    async def _cleanup_process(self, proc, communication_task, target_log):
        cleanup_task = asyncio.create_task(
            self._terminate_process_group(proc, communication_task, target_log)
        )
        try:
            await asyncio.shield(cleanup_task)
        except asyncio.CancelledError:
            await cleanup_task
            raise

    async def _start_process(self, script_path, target, environ):
        create_task = asyncio.create_task(
            asyncio.create_subprocess_exec(
                sys.executable,
                "-mgwf.exec",
                script_path,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=target.working_dir,
                env=environ,
                start_new_session=True,
            )
        )
        try:
            return await asyncio.shield(create_task)
        except asyncio.CancelledError:
            try:
                proc = await create_task
            except Exception:
                pass
            else:
                communication_task = asyncio.create_task(proc.communicate())
                await self._cleanup_process(proc, communication_task, target.log)
            raise

    async def _execute_task(self, target):
        with tempfile.NamedTemporaryFile(mode="w", prefix="gwf_") as script_file:
            script_file.write(target.script)
            script_file.flush()
            environ = os.environ.copy()
            environ["GWF_EXEC_WORKFLOW_ROOT"] = str(self.working_dir)
            proc = await self._start_process(
                script_file.name,
                target,
                environ,
            )
            target.log.debug(
                "started process %d in %s using script %s",
                proc.pid,
                target.working_dir,
                script_file.name,
            )
            communication_task = asyncio.create_task(proc.communicate())
            try:
                if target.time_limit is None:
                    stdout, stderr = await asyncio.shield(communication_task)
                else:
                    stdout, stderr = await asyncio.wait_for(
                        asyncio.shield(communication_task),
                        timeout=target.time_limit,
                    )
            except asyncio.TimeoutError:
                await self._cleanup_process(proc, communication_task, target.log)
                raise TimeLimitExceededError()
            except asyncio.CancelledError:
                await self._cleanup_process(proc, communication_task, target.log)
                raise

        # TODO: This should be made streaming..
        stdout_path, stderr_path = get_log_paths(self.working_dir, target.name)
        target.log.debug(
            "writing %d stdout byte(s) to %s and %d stderr byte(s) to %s",
            len(stdout),
            stdout_path,
            len(stderr),
            stderr_path,
        )
        with open(stdout_path, "wb") as log_file:
            log_file.write(stdout)
        with open(stderr_path, "wb") as log_file:
            log_file.write(stderr)

        if proc.returncode is not None and proc.returncode != 0:
            raise TaskFailedError(proc.returncode)


@attrs.define
class Server:
    scheduler: Scheduler = attrs.field()
    connections: set = attrs.field(factory=set)
    stop_requested: asyncio.Event = attrs.field(factory=asyncio.Event)

    async def send_response(self, writer, kind, **kwargs):
        writer.write(encode(kind, **kwargs).encode("utf-8"))
        await writer.drain()

    async def handle_message(self, kind, message):
        close_connection = False
        stop_server = False

        if kind == "enqueue_task":
            task = {
                "name": message.pop("name"),
                "script": message.pop("script"),
                "time_limit": message.pop("time_limit", None),
                "working_dir": message.pop("working_dir"),
                "deps": message.pop("deps"),
                "cores": message.pop("cores", 1),
                "memory": message.pop("memory", "1g"),
            }
            try:
                tid = self.scheduler.enqueue_task(**task)
            except (BackendError, ValueError) as error:
                response = ("task_rejected", {"message": str(error)})
            else:
                response = ("task_enqueued", {"tid": tid})
        elif kind == "get_task_states":
            response = (
                "task_states",
                {"tasks": self.scheduler.get_task_states()},
            )
        elif kind == "cancel_task":
            tid = message.pop("tid")
            await self.scheduler.cancel_task(tid)
            response = ("task_cancelled", {"tid": tid})
        elif kind == "shutdown":
            await self.scheduler.shutdown()
            response = ("shutdown", {})
            close_connection = True
            stop_server = True
        elif kind == "close":
            response = None
            close_connection = True
        else:
            raise ValueError(f"Unknown message kind: {kind}")

        assert not message, f"message of kind {kind} has not been fully parsed"
        return response, close_connection, stop_server

    async def handle_connection(self, reader, writer):
        stop_server = False
        self.connections.add(writer)
        try:
            while True:
                data = await reader.readline()
                if not data:
                    break
                message = json.loads(data)

                kind = message.pop("__kind__")
                response, close_connection, should_stop = await self.handle_message(
                    kind,
                    message,
                )
                stop_server = stop_server or should_stop
                if response is not None:
                    response_kind, response_data = response
                    await self.send_response(
                        writer,
                        response_kind,
                        **response_data,
                    )
                if close_connection:
                    break
        finally:
            self.connections.discard(writer)
            writer.close()
            with suppress(ConnectionError):
                await writer.wait_closed()
            if stop_server:
                self.stop_requested.set()

    async def start_server(self, socket_path, ready_event=None):
        socket_path = Path(socket_path)
        with _socket_lock(socket_path):
            socket_path.unlink(missing_ok=True)
            server = await asyncio.start_unix_server(
                self.handle_connection,
                path=socket_path,
            )
            socket_path.chmod(0o600)
            logger.info(
                "Started local scheduler with %d core(s) and %s memory, "
                "listening at %s",
                self.scheduler.max_cores,
                format_memory(self.scheduler.max_memory),
                socket_path,
            )
            if ready_event is not None:
                ready_event.set()
            try:
                async with server:
                    try:
                        await self.stop_requested.wait()
                    except asyncio.CancelledError:
                        logger.info("Shutting down...")
                    finally:
                        connections = list(self.connections)
                        for writer in connections:
                            writer.close()
                        await asyncio.gather(
                            *(writer.wait_closed() for writer in connections),
                            return_exceptions=True,
                        )
            finally:
                if not self.stop_requested.is_set():
                    await self.scheduler.shutdown()
                socket_path.unlink(missing_ok=True)
                logger.info("Bye!")


async def start_cluster_async(
    working_dir,
    max_cores,
    max_memory=None,
    ready_event=None,
):
    scheduler = Scheduler(working_dir, max_cores, max_memory)
    s = Server(scheduler)
    await s.start_server(get_socket_path(working_dir), ready_event)


def start_cluster(*args, debug=False, **kwargs):
    asyncio.run(start_cluster_async(*args, **kwargs), debug=debug)


setup = (create_backend, 0)
