import pytest

from gwf.backends.base import BackendError
from gwf.backends.logmanager import MemoryLogManager
from gwf.backends.local import (
    Client,
    Executor,
    LocalBackend,
    LocalStatus,
    TaskScheduler,
    Task,
)
from gwf import Target
from gwf.backends import Status
from gwf.backends.exceptions import DependencyError


class Fakescheduler:
    def __init__(self):
        self.history = []
        self._task_status = {}

    def enqueue_task(self, task):
        self.set_status(task.id, LocalStatus.SUBMITTED)

    def cancel_task(self, task_id):
        self.set_status(task_id, LocalStatus.CANCELLED)

    def set_status(self, task_id, status):
        self._task_status[task_id] = status
        self.history.append((task_id, status))

    def get_status(self, task_id):
        return self._task_status[task_id]


@pytest.fixture
def scheduler():
    return Fakescheduler()


class TestExecutor:
    def test_task_successful(self, scheduler, log_manager):
        executor = Executor(scheduler, log_manager=log_manager)
        executor.execute(Task(id="foo", script="exit 0"))
        executor.wait()

        assert scheduler.history[0] == ("foo", LocalStatus.RUNNING)
        assert scheduler.history[1] == ("foo", LocalStatus.COMPLETED)

    def test_task_failed(self, scheduler, log_manager):
        executor = Executor(scheduler, log_manager=log_manager)
        executor.execute(Task(id="foo", script="exit 1"))
        executor.wait()

        assert scheduler.history[0] == ("foo", LocalStatus.RUNNING)
        assert scheduler.history[1] == ("foo", LocalStatus.FAILED)

    def test_cancel(self, scheduler, log_manager):
        executor = Executor(scheduler, log_manager=log_manager)
        executor.execute(Task(id="foo", script="sleep 1"))
        executor.cancel()
        executor.wait()

        assert scheduler.history[0] == ("foo", LocalStatus.RUNNING)
        assert scheduler.history[1] == ("foo", LocalStatus.CANCELLED)

    def test_terminate(self, scheduler, log_manager):
        executor = Executor(scheduler, log_manager=log_manager)
        executor.execute(Task(id="foo", script="sleep 1"))
        executor.terminate()
        executor.wait()

        assert scheduler.history[0] == ("foo", LocalStatus.RUNNING)
        assert scheduler.history[1] == ("foo", LocalStatus.FAILED)


class TestScheduler:
    def test_task_lifecycle_successful(self, log_manager):
        scheduler = TaskScheduler(max_cores=1, log_manager=MemoryLogManager())

        task = Task(id="foo", script="sleep 1")
        scheduler.enqueue_task(task)
        assert scheduler.get_status("foo") == LocalStatus.SUBMITTED

        scheduler.schedule_once()
        assert scheduler.get_status("foo") == LocalStatus.RUNNING

        scheduler.wait()
        assert scheduler.get_status("foo") == LocalStatus.COMPLETED

    def test_cancel_task(self, log_manager):
        scheduler = TaskScheduler(max_cores=1, log_manager=MemoryLogManager())

        task = Task(id="foo", script="sleep 10")
        scheduler.enqueue_task(task)
        assert scheduler.get_status("foo") == LocalStatus.SUBMITTED

        scheduler.schedule_once()
        assert scheduler.get_status("foo") == LocalStatus.RUNNING

        scheduler.cancel_task("foo")
        scheduler.wait()
        assert scheduler.get_status("foo") == LocalStatus.CANCELLED

    def test_cap_at_one_core(self, log_manager):
        scheduler = TaskScheduler(max_cores=1, log_manager=log_manager)

        task1 = Task(id="foo1", script="sleep 1")
        scheduler.enqueue_task(task1)

        task2 = Task(id="foo2", script="sleep 1")
        scheduler.enqueue_task(task2)

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.RUNNING
        assert scheduler.get_status("foo2") == LocalStatus.SUBMITTED

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.SUBMITTED

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.RUNNING

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED

    def test_cap_at_two_cores(self, log_manager):
        scheduler = TaskScheduler(max_cores=2, log_manager=log_manager)

        task1 = Task(id="foo1", script="sleep 1")
        scheduler.enqueue_task(task1)

        task2 = Task(id="foo2", script="sleep 1")
        scheduler.enqueue_task(task2)

        task3 = Task(id="foo3", script="sleep 1")
        scheduler.enqueue_task(task3)

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.RUNNING
        assert scheduler.get_status("foo2") == LocalStatus.RUNNING
        assert scheduler.get_status("foo3") == LocalStatus.SUBMITTED

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo3") == LocalStatus.SUBMITTED

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo3") == LocalStatus.RUNNING

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo3") == LocalStatus.COMPLETED

    def test_unknown_dependency(self, log_manager):
        scheduler = TaskScheduler(max_cores=2, log_manager=log_manager)

        task1 = Task(id="foo1", script="sleep 1", dependencies=set(["bar"]))
        with pytest.raises(BackendError):
            scheduler.enqueue_task(task1)

    def test_wait_for_dependency(self, log_manager):
        scheduler = TaskScheduler(max_cores=2, log_manager=log_manager)

        task1 = Task(id="foo1", script="sleep 1")
        scheduler.enqueue_task(task1)

        task2 = Task(id="foo2", script="sleep 1", dependencies=set(["foo1"]))
        scheduler.enqueue_task(task2)

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.RUNNING
        assert scheduler.get_status("foo2") == LocalStatus.SUBMITTED

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.SUBMITTED

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.RUNNING

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED

    def test_wait_for_multiple_dependencies(self, log_manager):
        scheduler = TaskScheduler(max_cores=3, log_manager=log_manager)

        task1 = Task(id="foo1", script="sleep 1")
        scheduler.enqueue_task(task1)

        task2 = Task(id="foo2", script="sleep 1")
        scheduler.enqueue_task(task2)

        task3 = Task(id="foo3", script="sleep 1", dependencies=set(["foo1", "foo2"]))
        scheduler.enqueue_task(task3)

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.RUNNING
        assert scheduler.get_status("foo2") == LocalStatus.RUNNING
        assert scheduler.get_status("foo3") == LocalStatus.SUBMITTED

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo3") == LocalStatus.SUBMITTED

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo3") == LocalStatus.RUNNING

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo3") == LocalStatus.COMPLETED

    def test_dependents_fail_when_task_fails(self, log_manager):
        scheduler = TaskScheduler(max_cores=2, log_manager=log_manager)

        task1 = Task(id="foo1", script="exit 1")
        scheduler.enqueue_task(task1)

        task2 = Task(id="foo2", script="sleep 1")
        scheduler.enqueue_task(task2)

        task3 = Task(id="foo3", script="sleep 1", dependencies=set(["foo1", "foo2"]))
        scheduler.enqueue_task(task3)

        scheduler.schedule_once()
        assert scheduler.get_status("foo1") == LocalStatus.RUNNING
        assert scheduler.get_status("foo2") == LocalStatus.RUNNING
        assert scheduler.get_status("foo3") == LocalStatus.SUBMITTED

        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.FAILED
        assert scheduler.get_status("foo2") == LocalStatus.COMPLETED
        assert scheduler.get_status("foo3") == LocalStatus.FAILED

    def test_enqueue_task_with_failed_dependency(self, log_manager):
        scheduler = TaskScheduler(max_cores=2, log_manager=log_manager)

        task1 = Task(id="foo1", script="exit 1")
        scheduler.enqueue_task(task1)
        scheduler.schedule_once()
        scheduler.wait()
        assert scheduler.get_status("foo1") == LocalStatus.FAILED

        task2 = Task(id="foo2", script="sleep 1", dependencies=set(["foo1"]))
        scheduler.enqueue_task(task2)
        scheduler.schedule_once()
        assert scheduler.get_status("foo2") == LocalStatus.FAILED


class TestServer:
    def test_connect_close(self, local_backend):
        client = Client()
        client.connect()
        client.close()

    def test_submit_task_request(self, local_backend):
        client = Client()
        client.connect()
        task_id = client.submit(script="exit 0", working_dir="/tmp")
        assert client.status(task_id) == LocalStatus.SUBMITTED
        client.close()

    def test_cancel_task_request(self, local_backend):
        client = Client()
        client.connect()
        task_id = client.submit(script="sleep 5", working_dir="/tmp")
        client.cancel(task_id)
        assert client.status(task_id) == LocalStatus.CANCELLED
        client.close()


class TestLocalBackend:
    def test_connect_to_server_fails(self):
        with pytest.raises(BackendError):
            LocalBackend()

    def test_submit(self, local_backend):
        target = Target(
            "TestTarget", inputs=[], outputs=[], options={}, working_dir="/tmp"
        )
        with LocalBackend() as backend:
            assert backend.status(target) == Status.UNKNOWN
            backend.submit(target, dependencies=set())
            assert backend.status(target) == Status.SUBMITTED

    def test_submit_with_unknown_dependency(self, local_backend):
        target1 = Target(
            "TestTarget1", inputs=[], outputs=[], options={}, working_dir="/tmp"
        )
        target2 = Target(
            "TestTarget2", inputs=[], outputs=[], options={}, working_dir="/tmp"
        )
        with LocalBackend() as backend:
            with pytest.raises(DependencyError):
                backend.submit(target1, dependencies=set([target2]))

    def test_cancel(self, local_backend):
        target = Target(
            "TestTarget", inputs=[], outputs=[], options={}, working_dir="/tmp"
        )
        with LocalBackend() as backend:
            backend.submit(target, dependencies=set())
            assert backend.status(target) == Status.SUBMITTED
            backend.cancel(target)
            assert backend.status(target) == Status.UNKNOWN
