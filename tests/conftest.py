import functools
import time

import pytest

import gwf.conf
from gwf.backends.base import Backend, Status
from gwf.core import Graph, Target
from gwf.core import schedule as _schedule


@pytest.fixture(autouse=True)
def no_version_check(request, monkeypatch):
    monkeypatch.setitem(gwf.conf.CONFIG_DEFAULTS, "check_updates", False)


@pytest.fixture
def no_sleep(request, monkeypatch):
    def sleep(seconds):
        pass

    monkeypatch.setattr(time, "sleep", sleep)


class FakeBackend(Backend):
    option_defaults = {
        "cores": 1,
        "memory": "1g",
    }

    def __init__(self):
        super().__init__()
        self._tracked = {}

    def submit(self, target, dependencies):
        self._tracked[target] = Status.SUBMITTED

    def cancel(self, target):
        del self._tracked[target]

    def status(self, target):
        return self._tracked.get(target, Status.UNKNOWN)

    def close(self):
        pass

    def set_status(self, target, status):
        self._tracked[target] = status


class FakeFilesystem:
    def __init__(self):
        self._files = {}

    def add_file(self, path, changed_at):
        self._files[path] = changed_at

    def exists(self, path):
        return path in self._files

    def changed_at(self, path):
        if path not in self._files:
            raise FileNotFoundError(path)
        return self._files[path]


@pytest.fixture
def backend():
    return FakeBackend()


@pytest.fixture
def filesystem():
    return FakeFilesystem()


@pytest.fixture
def schedule(filesystem):
    return functools.partial(_schedule, filesystem=filesystem)


@pytest.fixture
def graph_factory():
    def factory(targets):
        return Graph.from_targets(targets)

    return factory


@pytest.fixture
def empty_graph(graph_factory):
    return graph_factory([])


@pytest.fixture
def trivial_graph(graph_factory):
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )
    return graph_factory([target])


@pytest.fixture
def diamond_graph(graph_factory):
    target1 = Target(
        "TestTarget1",
        inputs=[],
        outputs=["test_output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    target2 = Target(
        "TestTarget2",
        inputs=["test_output1.txt"],
        outputs=["test_output2.txt"],
        options={},
        working_dir="/some/dir",
    )
    target3 = Target(
        "TestTarget3",
        inputs=["test_output1.txt"],
        outputs=["test_output3.txt"],
        options={},
        working_dir="/some/dir",
    )
    target4 = Target(
        "TestTarget4",
        inputs=["test_output2.txt", "test_output3.txt"],
        outputs=["final_output.txt"],
        options={},
        working_dir="/some/dir",
    )
    return graph_factory([target1, target2, target3, target4])
