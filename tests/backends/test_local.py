import asyncio

import pytest
import pytest_asyncio

from gwf.backends.local import Client, Scheduler


@pytest_asyncio.fixture
async def s(tmp_path):
    tmp_path.joinpath(".gwf").mkdir()
    tmp_path.joinpath(".gwf", "logs").mkdir()
    return Scheduler(working_dir=tmp_path, max_cores=1)


@pytest.mark.asyncio
async def test_successful_task_without_deps(s):
    tid = await s.enqueue_task("foo", "exit 0", ".", None, set())
    await s.wait_for({tid}, timeout=1)
    assert s.get_task_state(tid) == "completed"


@pytest.mark.asyncio
async def test_successful_task_with_dependent(s):
    tid1 = await s.enqueue_task("foo", "exit 0", ".", None, set())
    tid2 = await s.enqueue_task("foo", "exit 0", ".", None, set([tid1]))
    await s.wait_for({tid1, tid2}, timeout=1)
    assert s.get_task_state(tid1) == "completed"
    assert s.get_task_state(tid2) == "completed"


@pytest.mark.asyncio
async def test_task_with_dependent_submitted_later(s):
    tid1 = await s.enqueue_task("foo", "exit 0", ".", None, set())
    await s.wait_for({tid1})
    assert s.get_task_state(tid1) == "completed"

    tid2 = await s.enqueue_task("foo", "exit 0", ".", None, set([tid1]))
    await s.wait_for({tid2})
    assert s.get_task_state(tid2) == "completed"


@pytest.mark.asyncio
async def test_failing_task_without_deps(s):
    tid = await s.enqueue_task("foo", "exit 1", ".", None, set())
    await s.wait_for({tid}, timeout=1)
    assert s.get_task_state(tid) == "failed"


@pytest.mark.asyncio
async def test_failed_task_with_dependents_1(s):
    tid1 = await s.enqueue_task("foo", "exit 1", ".", None, set())
    tid2 = await s.enqueue_task("foo", "exit 0", ".", None, set([tid1]))
    tid3 = await s.enqueue_task("foo", "exit 0", ".", None, set([tid2]))
    await asyncio.sleep(0.1)
    await s.cancel_task(tid1)
    await s.wait_for({tid1, tid2, tid3})
    assert s.get_task_state(tid1) == "failed"
    assert s.get_task_state(tid2) == "failed"
    assert s.get_task_state(tid3) == "failed"


@pytest.mark.asyncio
async def test_failed_task_with_dependents_2(s):
    tid1 = await s.enqueue_task("foo", "exit 0", ".", None, set())
    tid2 = await s.enqueue_task("foo", "exit 1", ".", None, set([tid1]))
    tid3 = await s.enqueue_task("foo", "exit 0", ".", None, set([tid2]))
    await asyncio.sleep(0.1)
    await s.cancel_task(tid1)
    await s.wait_for({tid1, tid2, tid3})
    assert s.get_task_state(tid1) == "completed"
    assert s.get_task_state(tid2) == "failed"
    assert s.get_task_state(tid3) == "failed"


@pytest.mark.asyncio
async def test_task_without_deps_times_out(s):
    tid = await s.enqueue_task("foo", "sleep 2", ".", 1, set())
    await s.wait_for({tid})
    assert s.get_task_state(tid) == "killed"


@pytest.mark.asyncio
async def test_task_without_deps_completes_within_timelimit(s):
    tid = await s.enqueue_task("foo", "sleep 1", ".", 2, set())
    await s.wait_for({tid})
    assert s.get_task_state(tid) == "completed"


@pytest.mark.asyncio
async def test_cancelled_task_without_deps(s):
    tid = await s.enqueue_task("foo", "sleep 3", ".", None, set())
    await asyncio.sleep(0.1)
    await s.cancel_task(tid)
    await s.wait_for({tid}, timeout=1)
    assert s.get_task_state(tid) == "cancelled"


@pytest.mark.asyncio
async def test_cancelled_task_with_dependents(s):
    tid1 = await s.enqueue_task("foo", "sleep 3", ".", None, set())
    tid2 = await s.enqueue_task("foo", "sleep 3", ".", None, set([tid1]))
    tid3 = await s.enqueue_task("foo", "sleep 3", ".", None, set([tid2]))
    await asyncio.sleep(0.1)
    await s.cancel_task(tid1)
    await s.wait_for({tid1, tid2, tid3})
    assert s.get_task_state(tid1) == "cancelled"
    assert s.get_task_state(tid2) == "cancelled"
    assert s.get_task_state(tid3) == "cancelled"


@pytest.mark.asyncio
async def test_task_writes_log_file(s):
    tid = await s.enqueue_task("foo", "echo hello world", ".", None, set())
    await s.wait_for({tid})
    contents = s.working_dir.joinpath(".gwf", "logs", "foo.stdout").read_text()
    assert contents == "hello world\n"


def test_client_connection_failure():
    with pytest.raises(Exception):
        Client.connect("localhost", 54321, attempts=1)
