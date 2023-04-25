import threading
import time

import attrs
import pytest

from gwf.backends.base import Backend, BackendStatus
from gwf.backends.local import Cluster
from gwf.core import Graph, Target, hash_spec


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
        self._tracked[target] = BackendStatus.SUBMITTED

    def cancel(self, target):
        del self._tracked[target]

    def status(self, target):
        return self._tracked.get(target, BackendStatus.UNKNOWN)

    def close(self):
        pass

    def set_status(self, target, status):
        self._tracked[target] = status


@attrs.define
class FakeFilesystem:
    files: dict = attrs.field(factory=dict)

    def add_file(self, path, changed_at):
        self.files[path] = changed_at

    def exists(self, path):
        return path in self.files

    def changed_at(self, path):
        if path not in self.files:
            raise FileNotFoundError(path)
        return self.files[path]


@attrs.define
class FakeSpecHashes:
    hashes: dict = attrs.field(factory=dict)

    def has_changed(self, target):
        spec_hash = hash_spec(target.spec)
        saved_hash = self.hashes.get(target)
        if saved_hash is None or spec_hash != saved_hash:
            return spec_hash
        return None

    def update(self, target):
        self.hashes[target] = hash_spec(target.spec)

    def invalidate(self, target):
        del self.hashes[target]

    def clear(self):
        self.hashes.clear()

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


@pytest.fixture
def backend():
    return FakeBackend()


@pytest.fixture(scope="module")
def local_backend():
    cluster = Cluster(num_workers=1)
    thread = threading.Thread(target=cluster.start)
    thread.start()
    yield
    cluster.shutdown()
    thread.join()


@pytest.fixture
def filesystem():
    return FakeFilesystem()


@pytest.fixture
def spec_hashes():
    return FakeSpecHashes()


@pytest.fixture
def empty_graph(graph_factory):
    return graph_factory([])


@pytest.fixture
def trivial_graph(filesystem):
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )
    return Graph.from_targets([target], filesystem)


@pytest.fixture
def diamond_graph(filesystem):
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
    return Graph.from_targets([target1, target2, target3, target4], filesystem)


@pytest.fixture
def simple_workflow(tmpdir):
    workflow_file = tmpdir.join("workflow.py")
    workflow_file.write(
        """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt'])
gwf.target('Target2', inputs=['a.txt'], outputs=['b.txt'])
gwf.target('Target3', inputs=['a.txt'], outputs=['c.txt'])
"""
    )
    with tmpdir.as_cwd():
        yield tmpdir


@pytest.fixture
def linear_workflow(tmpdir):
    workflow_file = tmpdir.join("workflow.py")
    workflow_file.write(
        """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=['a.txt'], outputs=['b.txt'])
gwf.target('Target2', inputs=['b.txt'], outputs=['c.txt'])
gwf.target('Target3', inputs=['c.txt'], outputs=['d.txt'])
"""
    )
    with tmpdir.as_cwd():
        yield tmpdir
