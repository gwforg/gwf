import asyncio
import io
import logging
import os
from pathlib import Path
import shlex
import socket
import tempfile
import time

import pytest
import pytest_asyncio

from gwf.backends.exceptions import BackendError
from gwf.backends.local import (
    Client,
    LocalStatus,
    Scheduler,
    Server,
    decode,
    encode,
    get_socket_path,
)
from gwf.core import Target
from gwf.executors import serialize
from gwf.log_storage import init_log_storage, prepare_log_storage_for_target
from tests.local_backend import background_scheduler


@pytest.fixture
def runtime_directory(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="gwf-test-", dir="/tmp") as directory:
        monkeypatch.setenv("XDG_RUNTIME_DIR", directory)
        yield Path(directory)


@pytest_asyncio.fixture
async def s(tmp_path):
    tmp_path.joinpath(".gwf").mkdir()
    init_log_storage(tmp_path)
    prepare_log_storage_for_target(tmp_path, "foo")
    scheduler = Scheduler(working_dir=tmp_path, max_cores=1, max_memory="8g")
    try:
        yield scheduler
    finally:
        await asyncio.wait_for(scheduler.shutdown(), timeout=15)


async def enqueue(
    s,
    name,
    spec,
    working_dir,
    time_limit,
    deps,
    cores=1,
    memory="1g",
):
    prepare_log_storage_for_target(s.working_dir, name)
    target = Target(
        name,
        inputs=[],
        outputs=[],
        options={},
        working_dir=working_dir,
        spec=spec,
    )
    script = io.StringIO()
    serialize(target, script)
    return s.enqueue_task(
        name,
        script.getvalue(),
        working_dir,
        time_limit,
        deps,
        cores,
        memory,
    )


async def wait_for_pids(path, timeout=5):
    async def _wait():
        while True:
            try:
                pids = [int(value) for value in path.read_text().split()]
            except (FileNotFoundError, ValueError):
                pass
            else:
                if len(pids) == 2 and all(pid > 0 for pid in pids):
                    return pids
            await asyncio.sleep(0.01)

    return await asyncio.wait_for(_wait(), timeout)


async def wait_for_path(path, timeout=5):
    async def _wait():
        while not path.exists():
            await asyncio.sleep(0.01)

    await asyncio.wait_for(_wait(), timeout)


async def wait_for_tasks(scheduler, tids, timeout=15):
    _, pending = await scheduler.wait_for(set(tids), timeout=timeout)
    assert not pending


async def wait_for_task(scheduler, tid, timeout=15):
    await wait_for_tasks(scheduler, {tid}, timeout=timeout)


async def assert_process_exits(pid, timeout=5):
    async def _wait():
        while True:
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                return
            await asyncio.sleep(0.01)

    await asyncio.wait_for(_wait(), timeout)


def process_tree_spec(pid_path):
    pid_path = shlex.quote(str(pid_path))
    return f"""
sleep 30 &
echo "$$ $!" > {pid_path}
wait
"""


def blocking_spec(started_path, release_path):
    started_path = shlex.quote(str(started_path))
    release_path = shlex.quote(str(release_path))
    return f"""
touch {started_path}
while [ ! -e {release_path} ]; do
    sleep 0.01
done
"""


@pytest.mark.asyncio
async def test_successful_task_without_deps(s):
    tid = await enqueue(s, "foo", "exit 0", ".", None, set())
    await wait_for_task(s, tid)
    assert s.get_task_state(tid) == LocalStatus.COMPLETED


@pytest.mark.asyncio
async def test_successful_task_with_dependent(s):
    tid1 = await enqueue(s, "foo", "exit 0", ".", None, set())
    tid2 = await enqueue(s, "foo", "exit 0", ".", None, set([tid1]))
    await wait_for_tasks(s, {tid1, tid2})
    assert s.get_task_state(tid1) == LocalStatus.COMPLETED
    assert s.get_task_state(tid2) == LocalStatus.COMPLETED


@pytest.mark.asyncio
async def test_task_with_dependent_submitted_later(s):
    tid1 = await enqueue(s, "foo", "exit 0", ".", None, set())
    await wait_for_task(s, tid1)
    assert s.get_task_state(tid1) == LocalStatus.COMPLETED

    tid2 = await enqueue(s, "foo", "exit 0", ".", None, set([tid1]))
    await wait_for_task(s, tid2)
    assert s.get_task_state(tid2) == LocalStatus.COMPLETED


@pytest.mark.asyncio
async def test_failing_task_without_deps(s):
    tid = await enqueue(s, "foo", "exit 1", ".", None, set())
    await wait_for_task(s, tid)
    assert s.get_task_state(tid) == LocalStatus.FAILED


@pytest.mark.asyncio
async def test_subprocess_start_failure_releases_core(s, tmp_path):
    tid = await enqueue(
        s,
        "foo",
        "exit 0",
        str(tmp_path / "missing"),
        None,
        set(),
    )
    await wait_for_task(s, tid)
    assert s.get_task_state(tid) == LocalStatus.FAILED

    tid = await enqueue(s, "foo", "exit 0", str(tmp_path), None, set())
    await wait_for_task(s, tid)
    assert s.get_task_state(tid) == LocalStatus.COMPLETED


@pytest.mark.asyncio
async def test_scheduler_allocates_requested_cores(tmp_path):
    scheduler = Scheduler(working_dir=tmp_path, max_cores=3, max_memory="8g")
    started = {
        name: tmp_path / f"{name}.started" for name in ("first", "second", "third")
    }
    releases = {
        name: tmp_path / f"{name}.release" for name in ("first", "second", "third")
    }

    try:
        first = await enqueue(
            scheduler,
            "first",
            blocking_spec(started["first"], releases["first"]),
            str(tmp_path),
            None,
            set(),
            cores=2,
        )
        second = await enqueue(
            scheduler,
            "second",
            blocking_spec(started["second"], releases["second"]),
            str(tmp_path),
            None,
            set(),
            cores=2,
        )
        third = await enqueue(
            scheduler,
            "third",
            blocking_spec(started["third"], releases["third"]),
            str(tmp_path),
            None,
            set(),
            cores=1,
        )

        await wait_for_path(started["first"])
        await wait_for_path(started["third"])
        assert scheduler.get_task_state(first) == LocalStatus.RUNNING
        assert scheduler.get_task_state(second) == LocalStatus.SUBMITTED
        assert scheduler.get_task_state(third) == LocalStatus.RUNNING
        assert scheduler.available_cores == 0

        releases["third"].touch()
        await wait_for_task(scheduler, third)
        assert scheduler.get_task_state(second) == LocalStatus.SUBMITTED
        assert scheduler.available_cores == 1

        releases["first"].touch()
        await wait_for_path(started["second"])
        assert scheduler.get_task_state(second) == LocalStatus.RUNNING
        assert scheduler.available_cores == 1

        releases["second"].touch()
        await wait_for_task(scheduler, second)
        assert scheduler.available_cores == 3
    finally:
        for release in releases.values():
            release.touch()
        await scheduler.shutdown()


@pytest.mark.asyncio
async def test_scheduler_allocates_requested_memory(tmp_path):
    scheduler = Scheduler(working_dir=tmp_path, max_cores=3, max_memory="3g")
    started = {
        name: tmp_path / f"{name}.started" for name in ("first", "second", "third")
    }
    releases = {
        name: tmp_path / f"{name}.release" for name in ("first", "second", "third")
    }

    try:
        first = await enqueue(
            scheduler,
            "first",
            blocking_spec(started["first"], releases["first"]),
            str(tmp_path),
            None,
            set(),
            memory="2g",
        )
        second = await enqueue(
            scheduler,
            "second",
            blocking_spec(started["second"], releases["second"]),
            str(tmp_path),
            None,
            set(),
            memory="2g",
        )
        third = await enqueue(
            scheduler,
            "third",
            blocking_spec(started["third"], releases["third"]),
            str(tmp_path),
            None,
            set(),
            memory="1g",
        )

        await wait_for_path(started["first"])
        await wait_for_path(started["third"])
        assert scheduler.get_task_state(first) == LocalStatus.RUNNING
        assert scheduler.get_task_state(second) == LocalStatus.SUBMITTED
        assert scheduler.get_task_state(third) == LocalStatus.RUNNING
        assert scheduler.available_memory == 0

        releases["third"].touch()
        await wait_for_task(scheduler, third)
        assert scheduler.get_task_state(second) == LocalStatus.SUBMITTED

        releases["first"].touch()
        await wait_for_path(started["second"])
        assert scheduler.get_task_state(second) == LocalStatus.RUNNING

        releases["second"].touch()
        await wait_for_task(scheduler, second)
        assert scheduler.available_memory == scheduler.max_memory
    finally:
        for release in releases.values():
            release.touch()
        await scheduler.shutdown()


@pytest.mark.asyncio
async def test_scheduler_rejects_target_requesting_too_many_cores(tmp_path):
    scheduler = Scheduler(working_dir=tmp_path, max_cores=1, max_memory="8g")

    with pytest.raises(BackendError, match="requested 2 cores"):
        await enqueue(
            scheduler,
            "foo",
            "",
            str(tmp_path),
            None,
            set(),
            cores=2,
        )

    assert scheduler.get_task_states() == {}


@pytest.mark.asyncio
async def test_scheduler_rejects_target_requesting_too_much_memory(tmp_path):
    scheduler = Scheduler(working_dir=tmp_path, max_cores=1, max_memory="1g")

    with pytest.raises(BackendError, match="requested 2.0GiB"):
        await enqueue(
            scheduler,
            "foo",
            "",
            str(tmp_path),
            None,
            set(),
            memory="2g",
        )

    assert scheduler.get_task_states() == {}


@pytest.mark.asyncio
async def test_scheduler_rejects_tasks_after_shutdown(tmp_path):
    scheduler = Scheduler(working_dir=tmp_path, max_cores=1, max_memory="1g")
    await scheduler.shutdown()

    with pytest.raises(BackendError, match="shutting down"):
        await enqueue(
            scheduler,
            "foo",
            "",
            str(tmp_path),
            None,
            set(),
        )


@pytest.mark.asyncio
async def test_failed_task_with_dependents_1(s):
    tid1 = await enqueue(s, "foo", "exit 1", ".", None, set())
    tid2 = await enqueue(s, "foo", "exit 0", ".", None, set([tid1]))
    tid3 = await enqueue(s, "foo", "exit 0", ".", None, set([tid2]))
    await wait_for_task(s, tid1)
    await wait_for_tasks(s, {tid1, tid2, tid3})
    assert s.get_task_state(tid1) == LocalStatus.FAILED
    assert s.get_task_state(tid2) == LocalStatus.FAILED
    assert s.get_task_state(tid3) == LocalStatus.FAILED


@pytest.mark.asyncio
async def test_failed_task_with_dependents_2(s):
    tid1 = await enqueue(s, "foo", "exit 0", ".", None, set())
    tid2 = await enqueue(s, "foo", "exit 1", ".", None, set([tid1]))
    tid3 = await enqueue(s, "foo", "exit 0", ".", None, set([tid2]))
    await wait_for_task(s, tid1)
    await wait_for_tasks(s, {tid1, tid2, tid3})
    assert s.get_task_state(tid1) == LocalStatus.COMPLETED
    assert s.get_task_state(tid2) == LocalStatus.FAILED
    assert s.get_task_state(tid3) == LocalStatus.FAILED


@pytest.mark.asyncio
async def test_task_without_deps_times_out(s, tmp_path):
    pid_path = tmp_path / "pids"
    tid = await enqueue(
        s,
        "foo",
        process_tree_spec(pid_path),
        str(tmp_path),
        2,
        set(),
    )
    pids = await wait_for_pids(pid_path)
    await wait_for_task(s, tid)

    assert s.get_task_state(tid) == LocalStatus.KILLED
    for pid in pids:
        await assert_process_exits(pid)


@pytest.mark.asyncio
async def test_task_without_deps_completes_within_timelimit(s):
    tid = await enqueue(s, "foo", "sleep 1", ".", 5, set())
    await wait_for_task(s, tid)
    assert s.get_task_state(tid) == LocalStatus.COMPLETED


@pytest.mark.asyncio
async def test_cancelled_task_without_deps(s, tmp_path):
    pid_path = tmp_path / "pids"
    tid = await enqueue(
        s,
        "foo",
        process_tree_spec(pid_path),
        str(tmp_path),
        None,
        set(),
    )
    pids = await wait_for_pids(pid_path)
    await asyncio.gather(s.cancel_task(tid), s.cancel_task(tid))

    assert s.get_task_state(tid) == LocalStatus.CANCELLED
    for pid in pids:
        await assert_process_exits(pid)


@pytest.mark.asyncio
async def test_cancelled_submitted_task_does_not_release_core(s, tmp_path):
    running_pid_path = tmp_path / "running-pids"
    running_tid = await enqueue(
        s,
        "foo",
        process_tree_spec(running_pid_path),
        str(tmp_path),
        None,
        set(),
    )
    await wait_for_pids(running_pid_path)

    cancelled_output = tmp_path / "cancelled"
    cancelled_tid = await enqueue(
        s,
        "foo",
        f"touch {shlex.quote(str(cancelled_output))}",
        str(tmp_path),
        None,
        set(),
    )
    await s.cancel_task(cancelled_tid)
    assert s.get_task_state(cancelled_tid) == LocalStatus.CANCELLED
    assert not cancelled_output.exists()

    queued_output = tmp_path / "queued"
    queued_tid = await enqueue(
        s,
        "foo",
        f"touch {shlex.quote(str(queued_output))}",
        str(tmp_path),
        None,
        set(),
    )
    await asyncio.sleep(0.2)
    assert s.get_task_state(queued_tid) == LocalStatus.SUBMITTED
    assert not queued_output.exists()

    await s.cancel_task(running_tid)
    await wait_for_task(s, queued_tid)
    assert s.get_task_state(queued_tid) == LocalStatus.COMPLETED
    assert queued_output.exists()


@pytest.mark.asyncio
async def test_cancelled_task_with_dependents(s):
    tid1 = await enqueue(s, "foo", "sleep 3", ".", None, set())
    tid2 = await enqueue(s, "foo", "sleep 3", ".", None, set([tid1]))
    tid3 = await enqueue(s, "foo", "sleep 3", ".", None, set([tid2]))
    await asyncio.sleep(0.1)
    await s.cancel_task(tid1)
    await wait_for_tasks(s, {tid1, tid2, tid3})
    assert s.get_task_state(tid1) == LocalStatus.CANCELLED
    assert s.get_task_state(tid2) == LocalStatus.CANCELLED
    assert s.get_task_state(tid3) == LocalStatus.CANCELLED


@pytest.mark.asyncio
async def test_task_writes_log_file(s):
    tid = await enqueue(s, "foo", "echo hello world", ".", None, set())
    await wait_for_task(s, tid)
    contents = s.working_dir.joinpath(
        ".gwf", "logs", "2", "c", "foo.stdout"
    ).read_text()
    assert contents == "hello world\n"


@pytest.mark.asyncio
async def test_task_debug_logs_include_target_context(s, caplog):
    with caplog.at_level(logging.DEBUG, logger="gwf.backends.local"):
        tid = await enqueue(s, "foo", "echo hello", ".", None, set())
        await wait_for_task(s, tid)

    messages = [
        record.getMessage()
        for record in caplog.records
        if record.name == "gwf.backends.local"
    ]
    target_messages = [
        message for message in messages if message.startswith("Target foo (id 0):")
    ]

    assert len(target_messages) >= 5
    assert any(
        "queued with 1 core(s), 1.0GiB memory, time limit None" in msg
        for msg in messages
    )
    assert any("starting with 1 core(s) and 1.0GiB memory" in msg for msg in messages)
    assert any("started process" in msg for msg in messages)
    assert any("writing 6 stdout byte(s)" in msg for msg in messages)
    assert any("state changed from running to completed" in msg for msg in messages)


@pytest.mark.asyncio
async def test_task_uses_gwf_exec(s, tmp_path):
    (tmp_path / "input.txt").write_text("input data")
    target = Target(
        "foo",
        inputs=["input.txt"],
        outputs=["output.txt"],
        options={},
        working_dir=str(tmp_path),
        spec="cat input.txt > output.txt",
    )
    script = io.StringIO()
    serialize(target, script)

    tid = s.enqueue_task(
        "foo",
        script.getvalue(),
        str(tmp_path),
        None,
        set(),
    )
    await wait_for_task(s, tid)

    assert s.get_task_state(tid) == LocalStatus.COMPLETED
    assert (tmp_path / "output.txt").read_text() == "input data"


def test_socket_path_uses_xdg_runtime_directory(tmp_path, runtime_directory):
    workflow_root = tmp_path / "workflow"
    other_workflow_root = tmp_path / "other-workflow"
    workflow_root.mkdir()
    other_workflow_root.mkdir()

    socket_path = get_socket_path(workflow_root)

    assert socket_path.parent == runtime_directory / "gwf"
    assert socket_path == get_socket_path(workflow_root / ".." / "workflow")
    assert socket_path != get_socket_path(other_workflow_root)


def test_local_cluster_socket_lifecycle(tmp_path, runtime_directory):
    with background_scheduler(tmp_path, 1) as socket_path:
        with Client.connect(socket_path) as client:
            assert client.status() == {}
        assert socket_path.is_socket()

        duplicate = Server(Scheduler(tmp_path, 1))
        with pytest.raises(BackendError, match="already running"):
            asyncio.run(duplicate.start_server(socket_path))

    assert not socket_path.exists()


def test_local_cluster_replaces_stale_socket(tmp_path, runtime_directory):
    socket_path = get_socket_path(tmp_path)
    stale_socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    stale_socket.bind(os.fspath(socket_path))
    stale_socket.close()

    with background_scheduler(tmp_path, 1) as socket_path:
        with Client.connect(socket_path) as client:
            assert client.status() == {}

    assert not socket_path.exists()


def test_client_submits_requested_cores(tmp_path, runtime_directory):
    target = Target(
        "foo",
        inputs=[],
        outputs=[],
        options={"cores": 2},
        working_dir=str(tmp_path),
        spec="exit 0",
    )

    with background_scheduler(tmp_path, 1) as socket_path:
        with Client.connect(socket_path) as client:
            with pytest.raises(BackendError, match="requested 2 cores"):
                client.submit(target)


def test_client_submits_requested_memory(tmp_path, runtime_directory):
    target = Target(
        "foo",
        inputs=[],
        outputs=[],
        options={"cores": 1, "memory": "2g"},
        working_dir=str(tmp_path),
        spec="exit 0",
    )

    with background_scheduler(tmp_path, 1, max_memory="1g") as socket_path:
        with Client.connect(socket_path) as client:
            with pytest.raises(BackendError, match="requested 2.0GiB"):
                client.submit(target)


def test_client_submits_target_walltime(tmp_path, runtime_directory):
    tmp_path.joinpath(".gwf").mkdir()
    init_log_storage(tmp_path)
    prepare_log_storage_for_target(tmp_path, "foo")
    target = Target(
        "foo",
        inputs=[],
        outputs=[],
        options={"cores": 1, "memory": "1m", "walltime": "00:00:01"},
        working_dir=str(tmp_path),
        spec="sleep 30",
    )

    with background_scheduler(tmp_path, 1, max_memory="1g") as socket_path:
        with Client.connect(socket_path) as client:
            tid = client.submit(target)
            deadline = time.monotonic() + 5
            while client.status()[str(tid)] != LocalStatus.KILLED:
                assert time.monotonic() < deadline
                time.sleep(0.05)


@pytest.mark.asyncio
async def test_server_rejects_tasks_after_shutdown_starts(
    tmp_path,
    runtime_directory,
):
    socket_path = get_socket_path(tmp_path)

    class BlockingScheduler:
        max_cores = 1
        max_memory = 1024**3

        def __init__(self):
            self.shutdown_started = asyncio.Event()
            self.allow_shutdown = asyncio.Event()
            self.task_enqueued = False
            self.accepting_tasks = True

        def enqueue_task(self, **kwargs):
            if not self.accepting_tasks:
                raise BackendError("Local scheduler is shutting down.")
            self.task_enqueued = True
            return 0

        async def shutdown(self):
            self.accepting_tasks = False
            self.shutdown_started.set()
            await self.allow_shutdown.wait()

    scheduler = BlockingScheduler()
    server = Server(scheduler)
    server_task = asyncio.create_task(server.start_server(socket_path))
    writers = []
    try:
        for _ in range(100):
            try:
                shutdown_reader, shutdown_writer = await asyncio.open_unix_connection(
                    socket_path
                )
                break
            except FileNotFoundError:
                await asyncio.sleep(0.01)
        else:
            pytest.fail("local server did not start")

        task_reader, task_writer = await asyncio.open_unix_connection(socket_path)
        writers.extend((shutdown_writer, task_writer))

        shutdown_writer.write(encode("shutdown").encode())
        await shutdown_writer.drain()
        await scheduler.shutdown_started.wait()

        task_writer.write(
            encode(
                "enqueue_task",
                name="late-task",
                script="",
                time_limit=None,
                working_dir=str(tmp_path),
                deps=[],
            ).encode()
        )
        await task_writer.drain()
        kind, response = decode(await task_reader.readline())

        assert kind == "task_rejected"
        assert response["message"] == "Local scheduler is shutting down."
        assert not scheduler.task_enqueued

        scheduler.allow_shutdown.set()
        kind, response = decode(await shutdown_reader.readline())
        assert kind == "shutdown"
        assert response == {}
        await asyncio.wait_for(server_task, timeout=5)
    finally:
        scheduler.allow_shutdown.set()
        for writer in writers:
            writer.close()
            await writer.wait_closed()
        if not server_task.done():
            server_task.cancel()
        await asyncio.gather(server_task, return_exceptions=True)


def test_client_connection_failure(tmp_path, runtime_directory):
    with pytest.raises(FileNotFoundError):
        Client.connect(get_socket_path(tmp_path))
