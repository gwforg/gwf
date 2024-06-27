"""Backend that runs targets on a local cluster.

To use this backend you must activate the `local` backend and start a local
cluster (with one or more workers) that the backend can submit targets to. To
start a cluster with two workers run the command::

    gwf -b local workers -n 2

in the working directory of your project. The workflow file must be accessible
to *gwf*. Thus, if your workflow file is not called `workflow.py` or the
workflow object is not called `gwf`, you must specify this so that *gwf* can
locate the workflow::

    gwf -f myworkflow.py:wf1 -b local workers -n 2

If the local backend is your default backend you can of course omit the
``-b local`` option.

If the ``-n`` option is omitted, *gwf* will detect the number of cores available
and use all of them.

To run your workflow, open another terminal and then type::

    gwf -b local run

To stop the pool of workers press :kbd:`Control-c`.

**Backend options:**

* **local.host (str):** Set the host that the workers are running on
    (default: localhost).
* **local.port (int):** Set the port used to connect to the workers
    (default: 12345).

**Target options:**

None available.
"""

import asyncio
import itertools
import json
import logging
import multiprocessing
import socket
import time
from enum import Enum
from io import TextIOWrapper
from pathlib import Path
from typing import Generator

import attrs

from .base import BackendStatus, TrackingBackend
from .exceptions import BackendError

__all__ = ("Client", "Server")


DEFAULT_HOST = "localhost"
DEFAULT_PORT = 12345

logger = logging.getLogger(__name__)


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
    def connect(cls, hostname=DEFAULT_HOST, port=DEFAULT_PORT, attempts=20):
        for attempts_used in range(attempts):
            try:
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.connect((hostname, port))
                return cls.from_socket(sock)
            except OSError:
                retry_delay = 2**attempts_used
                logger.warning(
                    "Could not connect, trying again in %d second(s). "
                    "Did you start workers with `gwf workers`?",
                    retry_delay,
                )
                time.sleep(retry_delay)
        raise ConnectionRefusedError("Failed to connect after three failed attempts")

    def submit(self, target, deps=None):
        self.send(
            "enqueue_task",
            name=target.name,
            script=target.spec,
            time_limit=None,
            working_dir=target.working_dir,
            deps=deps or [],
        )
        kind, response = self.recv()
        assert kind == "task_enqueued"
        return response["tid"]

    def status(self):
        self.send("get_task_states")
        msg_type, response = self.recv()
        assert msg_type == "task_states", "invalid response received"
        return {k: LocalStatus[v] for k, v in response["tasks"].items()}

    def cancel(self, job_id):
        self.send("cancel_task", tid=job_id)

    def shutdown(self):
        self.send("shutdown")

    def close(self):
        self.send("close")
        self.sock.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


@attrs.define
class LocalOps:
    working_dir: str = attrs.field()
    host: str = attrs.field()
    port: int = attrs.field()
    target_defaults: dict = attrs.field()

    _client: Client = attrs.field(init=False, repr=False)

    @_client.default
    def _create_client(self):
        try:
            return Client.connect(self.host, self.port)
        except ConnectionRefusedError:
            raise BackendError(
                f"Local backend could not connect to workers at {self.host} port {self.port}. "
            )

    def get_job_states(self, tracked_jobs):
        return {
            int(k): STATUS_MAP[v]
            for k, v in self._client.status().items()
            if int(k) in tracked_jobs
        }

    def submit_target(self, target, dependency_ids):
        return self._client.submit(
            target,
            deps=dependency_ids,
        )

    def cancel_job(self, job_id):
        self._client.cancel(job_id)

    def close(self):
        self._client.close()


def create_backend(working_dir, host=DEFAULT_HOST, port=DEFAULT_PORT):
    return TrackingBackend(
        working_dir,
        name="local",
        ops=LocalOps(working_dir, host, port, target_defaults={}),
    )


@attrs.frozen
class TaskFailedError(Exception):
    returncode: int


@attrs.frozen
class TimeLimitExceededError(Exception):
    pass


@attrs.define
class Scheduler:
    working_dir: Path = attrs.field(converter=Path)
    max_cores: int = attrs.field(default=multiprocessing.cpu_count())

    tid_generator: Generator = attrs.field(factory=itertools.count)
    events: asyncio.Queue = attrs.field(factory=asyncio.Queue)
    task_states: dict = attrs.field(factory=dict)
    tasks: dict = attrs.field(factory=dict)

    # Ressources
    cores_ressource: asyncio.Semaphore = attrs.field()

    @cores_ressource.default
    def create_cores_ressource(self):
        return asyncio.Semaphore(self.max_cores)

    async def enqueue_task(self, name, script, working_dir, time_limit, deps):
        tid = next(self.tid_generator)
        worker_task = asyncio.create_task(
            self.try_handle_task(
                tid,
                name,
                script,
                working_dir,
                time_limit,
                deps,
            )
        )
        self.tasks[tid] = worker_task
        self.task_states[tid] = LocalStatus.SUBMITTED
        return tid

    async def cancel_task(self, tid):
        if self.task_states[tid] in (LocalStatus.SUBMITTED, LocalStatus.RUNNING):
            worker_task = self.tasks[tid]
            worker_task.cancel()
            self.task_states[tid] = LocalStatus.CANCELLED

    async def kill(self):
        for worker in self.tasks.values():
            worker.cancel()

    async def wait(self):
        await asyncio.wait(self.tasks.values())

    async def shutdown(self):
        await self.kill()
        await self.wait()

    def get_task_state(self, tid):
        return self.task_states.get(tid, None)

    def get_task_states(self):
        return dict(self.task_states)

    async def wait_for(self, tids, timeout=None):
        tasks = {self.tasks[tid] for tid in tids}
        await asyncio.wait(tasks, timeout=timeout)

    async def _gentle_kill(self, proc):
        if proc is None:
            return

        proc.kill()
        await asyncio.sleep(1)
        if proc.returncode is None:
            await asyncio.sleep(10)
            proc.terminate()
        await proc.wait()

    async def try_handle_task(self, tid, name, script, working_dir, time_limit, deps):
        proc = None
        try:
            if deps:
                await asyncio.wait(
                    {self.tasks[tid] for tid in deps},
                    return_when=asyncio.ALL_COMPLETED,
                )
                for dep_tid in deps:
                    if self.task_states[dep_tid] != LocalStatus.COMPLETED:
                        self.task_states[tid] = self.task_states[dep_tid]
                        return

            await self.cores_ressource.acquire()
            self.task_states[tid] = LocalStatus.RUNNING

            proc = await asyncio.create_subprocess_shell(
                script,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=working_dir,
            )
            try:
                logger.debug("task starting")
                stdout, stderr = await asyncio.wait_for(
                    proc.communicate(), timeout=time_limit
                )
                logger.debug("writing log files")
                # TODO: This should be made streaming..
                with open(
                    self.working_dir.joinpath(".gwf", "logs", f"{name}.stdout"), "wb"
                ) as log_file:
                    log_file.write(stdout)
                with open(
                    self.working_dir.joinpath(".gwf", "logs", f"{name}.stderr"), "wb"
                ) as log_file:
                    log_file.write(stderr)
                logger.debug("wrote log files")
            except asyncio.TimeoutError:
                raise TimeLimitExceededError()

            if proc.returncode is not None and proc.returncode != 0:
                raise TaskFailedError(proc.returncode)
        except asyncio.CancelledError:
            logger.debug("got cancel for task")
            await self._gentle_kill(proc)
            self.task_states[tid] = LocalStatus.CANCELLED
            logger.debug("task states after cancel: %s", self.task_states)
        except TimeLimitExceededError:
            await self._gentle_kill(proc)
            self.task_states[tid] = LocalStatus.KILLED
        except TaskFailedError:
            self.task_states[tid] = LocalStatus.FAILED
        else:
            self.task_states[tid] = LocalStatus.COMPLETED
        finally:
            self.cores_ressource.release()


@attrs.define
class Server:
    scheduler: Scheduler = attrs.field()
    server: asyncio.base_events.Server = attrs.field(default=None)

    async def send_response(self, writer, kind, **kwargs):
        writer.write(encode(kind, **kwargs).encode("utf-8"))
        await writer.drain()

    async def handle_connection(self, reader, writer):
        while True:
            data = await reader.readline()
            if data is None:
                break
            message = json.loads(data)

            kind = message.pop("__kind__")
            if kind == "enqueue_task":
                tid = await self.scheduler.enqueue_task(
                    name=message.pop("name"),
                    script=message.pop("script"),
                    time_limit=message.pop("time_limit", None),
                    working_dir=message.pop("working_dir"),
                    deps=message.pop("deps"),
                )
                await self.send_response(
                    writer,
                    "task_enqueued",
                    tid=tid,
                )
            elif kind == "get_task_state":
                tid = message.pop("tid")
                await self.send_response(
                    writer, "task_state", state=self.scheduler.get_task_state(tid)
                )
            elif kind == "get_task_states":
                await self.send_response(
                    writer, "task_states", tasks=self.scheduler.get_task_states()
                )
            elif kind == "cancel_task":
                tid = message.pop("tid")
                await self.scheduler.cancel_task(tid)
            elif kind == "shutdown":
                self.server.close()
                await self.server.wait_closed()
                break
            elif kind == "close":
                break
            assert not message, f"message of kind {kind} has not been fully parsed"

    async def start_server(self, host=DEFAULT_HOST, port=DEFAULT_PORT):
        self.server = await asyncio.start_server(
            self.handle_connection,
            host,
            port,
            start_serving=False,
        )
        try:
            await self.server.serve_forever()
            addr, port = self.server.sockets[0].getsockname()
            logger.info("Listening on %s port %s", addr, port)
        except asyncio.CancelledError:
            logger.info("Shutting down...")
        finally:
            logger.info("Bye!")


async def start_cluster_async(working_dir, max_cores, host, port):
    scheduler = Scheduler(working_dir, max_cores)
    s = Server(scheduler)
    await s.start_server(host, port)


def start_cluster(*args, debug=False, **kwargs):
    asyncio.run(start_cluster_async(*args, **kwargs), debug=debug)


@attrs.frozen
class BackgroundCluster:
    host: str = attrs.field()
    port: int = attrs.field()
    process: multiprocessing.Process = attrs.field()

    @classmethod
    def start(cls, working_dir, max_cores, host, port, **kwargs):
        proc = multiprocessing.Process(
            target=start_cluster,
            args=(working_dir, max_cores, host, port),
            kwargs=kwargs,
        )
        proc.start()
        return cls(host, port, proc)

    def shutdown(self):
        with Client.connect("localhost", 12345) as c:
            c.shutdown()
        self.process.join(timeout=1)
        self.process.kill()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.shutdown()


def start_cluster_in_background(*args, **kwargs):
    return BackgroundCluster.start(*args, **kwargs)


setup = (create_backend, 0)
