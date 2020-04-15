import logging
import unittest
from unittest.mock import Mock, call, patch

import pytest

from gwf import AnonymousTarget, Graph, Scheduler, Target, Workflow
from gwf.backends import Backend, Status
from gwf.backends.exceptions import LogError
from gwf.exceptions import NameError, WorkflowError
from gwf.workflow import _flatten

FLATTEN_TESTS = [
    (["a", "b"], ["a", "b"]),
    ({"A": ["a1", "a2"]}, ["a1", "a2"]),
    ({"A": ["a1", "a2"], "B": ["b1", "b2"]}, ["a1", "a2", "b1", "b2"]),
    ([{"A": "a1", "B": "b1"}, {"A": "a2", "B": "b2"}], ["a1", "a2", "b1", "b2"]),
    (
        [{"A": ["a1", "A1"], "B": "b1"}, {"A": ["a2", "A2"], "B": "b2"}],
        ["a1", "a2", "b1", "b2", "A1", "A2"],
    ),
]


@pytest.mark.parametrize("value,expected", FLATTEN_TESTS)
def test_flatten(value, expected):
    assert set(_flatten(value)) == set(expected)


class DummyBackend(Backend):
    def __init__(self):
        super().__init__()
        self._tracked = {}

    def submit(self, target, dependencies):
        self._tracked[target] = Status.SUBMITTED

    def cancel(self, target):
        del self._tracked[target]

    def status(self, target):
        return self._tracked.get(target, Status.UNKNOWN)

    def logs(self, target, stderr=False):
        raise LogError

    def close(self):
        pass

    def set_status(self, target, status):
        assert status in (Status.RUNNING, Status.UNKNOWN)
        self._tracked[target] = status


@pytest.fixture
def backend():
    backend = DummyBackend()
    backend.submit = Mock(wraps=backend.submit)
    backend.status = Mock(wraps=backend.status)
    return backend


class TestTarget(unittest.TestCase):
    def test_target_with_invalid_name_raises_exception(self):
        with self.assertRaises(NameError):
            Target(
                "123abc", inputs=[], outputs=[], options={}, working_dir="/some/path"
            )

    def test_target_without_outputs_is_a_sink(self):
        target = Target(
            name="TestTarget",
            inputs=[],
            outputs=[],
            options={},
            working_dir="/some/path",
        )
        self.assertTrue(target.is_sink)

    def test_target_with_outputs_is_not_a_sink(self):
        target = Target(
            name="TestTarget",
            inputs=[],
            outputs=["test_output1.txt", "test_output2.txt"],
            options={},
            working_dir="/some/path",
        )
        self.assertFalse(target.is_sink)

    def test_target_without_inputs_is_a_source(self):
        target = Target(
            name="TestTarget",
            inputs=[],
            outputs=[],
            options={},
            working_dir="/some/path",
        )
        self.assertTrue(target.is_source)

    def test_target_with_inputs_is_not_a_source(self):
        target = Target(
            name="TestTarget",
            inputs=["test_input1.txt", "test_input2.txt"],
            outputs=[],
            options={},
            working_dir="/some/path",
        )
        self.assertFalse(target.is_source)

    def test_assigning_spec_to_target_sets_spec_attribute(self):
        target = (
            Target(
                name="TestTarget",
                inputs=[],
                outputs=[],
                options={},
                working_dir="/some/path",
            )
            << "this is a spec"
        )
        self.assertIsNotNone(target.spec)
        self.assertEqual(target.spec, "this is a spec")

    def test_inherit_options(self):
        target = Target(
            "TestTarget",
            inputs=[],
            outputs=[],
            options={"cores": 8},
            working_dir="/some/dir",
        )
        target.inherit_options({"cores": 4, "memory": "4g"})
        self.assertEqual(target.options, {"cores": 8, "memory": "4g"})

    def test_str_on_target(self):
        target = Target(
            "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/path"
        )
        self.assertEqual(str(target), "TestTarget")

    def test_input_is_empty_string(self):
        with self.assertRaises(WorkflowError):
            Target(
                name="TestTarget",
                inputs=["ab", ""],
                outputs=[],
                options={},
                working_dir="/some/path",
            )

    def test_output_is_empty_string(self):
        with self.assertRaises(WorkflowError):
            Target(
                name="TestTarget",
                inputs=[],
                outputs=["ab", ""],
                options={},
                working_dir="/some/path",
            )

    def test_input_contains_nonprintable(self):
        with self.assertRaises(WorkflowError):
            Target(
                name="TestTarget",
                inputs=["a\nb", "ac"],
                outputs=[],
                options={},
                working_dir="/some/path",
            )

    def test_output_contains_nonprintable(self):
        with self.assertRaises(WorkflowError):
            Target(
                name="TestTarget",
                inputs=[],
                outputs=["a\nb", "ac"],
                options={},
                working_dir="/some/path",
            )


def test_build_branch_join_graph():
    t1 = Target(
        name="Target1",
        inputs=["test_input1.txt", "test_input2.txt"],
        outputs=["t1_output1.txt", "t1_output2.txt"],
        options={},
        working_dir="/some/dir",
    )

    t2 = Target(
        name="Target2",
        inputs=["t1_output1.txt"],
        outputs=["t2_output.txt"],
        options={},
        working_dir="/some/dir",
    )

    t3 = Target(
        name="Target3",
        inputs=["t1_output2.txt"],
        outputs=["t3_output.txt"],
        options={},
        working_dir="/some/dir",
    )

    t4 = Target(
        name="Target4",
        inputs=["t2_output.txt", "t3_output.txt"],
        outputs=["t4_output.txt"],
        options={},
        working_dir="/some/dir",
    )

    targets = {"Target1": t1, "Target2": t2, "Target3": t3, "Target4": t4}

    graph = Graph.from_targets(targets)

    assert len(graph.targets) == 4

    assert not graph.dependencies[t1]
    assert graph.dependents[t1] == {t2, t3}

    assert graph.dependencies[t2] == {t1}
    assert graph.dependents[t2] == {t4}

    assert graph.dependencies[t3] == {t1}
    assert graph.dependents[t3] == {t4}

    assert graph.dependencies[t4] == {t2, t3}
    assert graph.dependents[t4] == set()

    assert graph.provides["/some/dir/t1_output1.txt"] == t1
    assert graph.provides["/some/dir/t1_output2.txt"] == t1
    assert graph.provides["/some/dir/t2_output.txt"] == t2
    assert graph.provides["/some/dir/t3_output.txt"] == t3
    assert graph.provides["/some/dir/t4_output.txt"] == t4

    assert graph.unresolved == {
        "/some/dir/test_input1.txt",
        "/some/dir/test_input2.txt",
    }


def test_graph_raises_multiple_providers_error():
    t1 = Target(
        name="Target1",
        inputs=[],
        outputs=["t1_output1.txt", "t1_output2.txt"],
        options={},
        working_dir="/some/dir",
    )

    t2 = Target(
        name="Target2",
        inputs=[],
        outputs=["t1_output2.txt", "t1_output3.txt"],
        options={},
        working_dir="/some/dir",
    )

    with pytest.raises(WorkflowError):
        Graph.from_targets({"Target1": t1, "Target2": t2})


def test_graph_raises_circular_dependency_error():
    t1 = Target(
        name="Target1",
        inputs=["f1.txt"],
        outputs=["f2.txt"],
        options={},
        working_dir="/some/dir",
    )
    t2 = Target(
        name="Target2",
        inputs=["f2.txt"],
        outputs=["f3.txt"],
        options={},
        working_dir="/some/dir",
    )
    t3 = Target(
        name="Target3",
        inputs=["f3.txt"],
        outputs=["f1.txt"],
        options={},
        working_dir="/some/dir",
    )
    with pytest.raises(WorkflowError):
        Graph.from_targets({"Target1": t1, "Target2": t2, "Target3": t3})


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


class TestShouldRun(unittest.TestCase):
    def setUp(self):
        workflow = Workflow(working_dir="/some/dir")
        self.target1 = workflow.target(
            "TestTarget1", inputs=[], outputs=["test_output1.txt"]
        )
        self.target2 = workflow.target(
            "TestTarget2", inputs=["test_output1.txt"], outputs=["test_output2.txt"]
        )
        self.target3 = workflow.target(
            "TestTarget3", inputs=["test_output1.txt"], outputs=["test_output3.txt"]
        )
        self.target4 = workflow.target(
            "TestTarget4",
            inputs=["test_output2.txt", "test_output3.txt"],
            outputs=["final_output.txt"],
        )

        self.graph = Graph.from_targets(workflow.targets)
        self.backend = DummyBackend()
        self.filesystem = FakeFilesystem()
        self.scheduler = Scheduler(
            graph=self.graph, backend=self.backend, filesystem=self.filesystem
        )

    def test_target_should_run_if_one_of_its_dependencies_does_not_exist(self):
        with self.assertLogs(level="DEBUG") as logs:
            self.assertTrue(self.scheduler.should_run(self.target1))

        self.assertEqual(
            logs.output,
            [
                "DEBUG:gwf.core:TestTarget1 should run because its output file /some/dir/test_output1.txt does not exist"
            ],
        )

    def test_target_should_run_if_one_of_its_dependencies_should_run(self):
        with self.assertLogs(level="DEBUG") as logs:
            self.assertTrue(self.scheduler.should_run(self.target2))

        self.assertEqual(
            logs.output,
            [
                "DEBUG:gwf.core:TestTarget1 should run because its output file /some/dir/test_output1.txt does not exist",
                "DEBUG:gwf.core:TestTarget2 should run because its dependency TestTarget1 should run",
            ],
        )

    def test_target_should_run_if_it_is_a_sink(self):
        target = Target(
            "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
        )
        graph = Graph.from_targets({"TestTarget": target})
        scheduler = Scheduler(
            graph=graph, backend=DummyBackend(), filesystem=self.filesystem
        )
        with self.assertLogs(level="DEBUG") as logs:
            self.assertTrue(scheduler.schedule(target))
            self.assertEqual(
                logs.output,
                [
                    "DEBUG:gwf.core:Scheduling target TestTarget",
                    "DEBUG:gwf.core:TestTarget should run because it is a sink",
                    "INFO:gwf.core:Submitting target TestTarget",
                ],
            )

    def test_target_should_not_run_if_it_is_a_source_and_all_outputs_exist(self):
        workflow = Workflow(working_dir="/some/dir")
        target = workflow.target(
            "TestTarget1", inputs=[], outputs=["test_output1.txt", "test_output2.txt"]
        )

        graph = Graph.from_targets(workflow.targets)

        self.filesystem.add_file("/some/dir/test_output1.txt", changed_at=1)
        self.filesystem.add_file("/some/dir/test_output2.txt", changed_at=2)

        scheduler = Scheduler(
            graph=graph, backend=DummyBackend(), filesystem=self.filesystem
        )
        self.assertFalse(scheduler.should_run(target))

    def test_should_run_if_any_input_file_is_newer_than_any_output_file(self):
        self.filesystem.add_file("/some/dir/test_output1.txt", changed_at=0)
        self.filesystem.add_file("/some/dir/test_output2.txt", changed_at=1)
        self.filesystem.add_file("/some/dir/test_output3.txt", changed_at=3)
        self.filesystem.add_file("/some/dir/final_output.txt", changed_at=2)

        self.assertFalse(self.scheduler.should_run(self.target1))
        self.assertFalse(self.scheduler.should_run(self.target2))
        self.assertFalse(self.scheduler.should_run(self.target3))
        self.assertTrue(self.scheduler.should_run(self.target4))

    def test_should_run_not_run_if_all_outputs_are_newer_then_the_inputs(self):
        self.filesystem.add_file("/some/dir/test_output1.txt", changed_at=0)
        self.filesystem.add_file("/some/dir/test_output2.txt", changed_at=1)
        self.filesystem.add_file("/some/dir/test_output3.txt", changed_at=3)
        self.filesystem.add_file("/some/dir/final_output.txt", changed_at=4)

        self.assertFalse(self.scheduler.should_run(self.target1))
        self.assertFalse(self.scheduler.should_run(self.target2))
        self.assertFalse(self.scheduler.should_run(self.target3))
        self.assertFalse(self.scheduler.should_run(self.target4))


def test_exception_if_input_file_is_not_provided_and_output_file_exists():
    workflow = Workflow(working_dir="/some/dir")
    target = workflow.target("TestTarget", inputs=["in.txt"], outputs=["out.txt"])
    graph = Graph.from_targets(workflow.targets)
    backend = DummyBackend()

    filesystem = FakeFilesystem()
    filesystem.add_file("/some/dir/out.txt", changed_at=1)

    scheduler = Scheduler(graph=graph, backend=backend, filesystem=filesystem)

    with pytest.raises(WorkflowError):
        scheduler.should_run(target)


@patch("gwf.core.os.path.exists", return_value=True, autospec=True)
def test_two_targets_producing_the_same_file_but_declared_with_rel_and_abs_path(
    mock_os_path_exists,
):
    workflow = Workflow(working_dir="/some/dir")
    workflow.target("TestTarget1", inputs=[], outputs=["/some/dir/test_output.txt"])
    workflow.target("TestTarget2", inputs=[], outputs=["test_output.txt"])

    with pytest.raises(WorkflowError):
        Graph.from_targets(workflow.targets)


def test_scheduling_submitted_target(backend, monkeypatch):
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )
    graph = Graph.from_targets({"TestTarget": target})
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    monkeypatch.setattr(scheduler, "should_run", lambda t: True)
    backend.submit(target, dependencies=set())
    assert len(backend.submit.call_args_list) == 1
    assert scheduler.schedule(target)
    assert len(backend.submit.call_args_list) == 1


def test_scheduling_unsubmitted_target(backend, monkeypatch):
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )
    graph = Graph.from_targets({"TestTarget": target})
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    monkeypatch.setattr(scheduler, "should_run", lambda t: True)
    assert scheduler.schedule(target)
    assert len(backend.submit.call_args_list) == 1
    assert call(target, dependencies=set()) in backend.submit.call_args_list


def test_non_existing_files_not_provided_by_other_target(backend):
    target = Target(
        "TestTarget",
        inputs=["test_input.txt"],
        outputs=[],
        options={},
        working_dir="/some/dir",
    )
    graph = Graph.from_targets({"TestTarget": target})
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    with pytest.raises(WorkflowError):
        scheduler.schedule(target)


def test_existing_files_not_provided_by_other_target(backend):
    target = Target(
        "TestTarget",
        inputs=["test_input.txt"],
        outputs=[],
        options={},
        working_dir="/some/dir",
    )
    graph = Graph.from_targets({"TestTarget": target})

    filesystem = FakeFilesystem()
    filesystem.add_file("/some/dir/test_input.txt", changed_at=0)

    scheduler = Scheduler(graph=graph, backend=backend, filesystem=filesystem)
    assert scheduler.schedule(target)


def test_scheduling_target_with_deps_that_are_not_submitted(backend, monkeypatch):
    target1 = Target(
        "TestTarget1",
        inputs=[],
        outputs=["test_output.txt"],
        options={},
        working_dir="/some/dir",
    )
    target2 = Target(
        "TestTarget2",
        inputs=["test_output.txt"],
        outputs=[],
        options={},
        working_dir="/some/dir",
    )
    graph = Graph.from_targets({"TestTarget1": target1, "TestTarget2": target2})
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    monkeypatch.setattr(scheduler, "should_run", lambda t: True)
    assert scheduler.schedule(target2)
    assert len(backend.submit.call_args_list) == 2
    assert call(target1, dependencies=set()) in backend.submit.call_args_list
    assert call(target2, dependencies=set([target1])) in backend.submit.call_args_list


def test_scheduling_target_with_deep_deps_that_are_not_submitted(backend, monkeypatch):
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
        inputs=["test_output2.txt"],
        outputs=["test_output3.txt"],
        options={},
        working_dir="/some/dir",
    )
    target4 = Target(
        "TestTarget4",
        inputs=["test_output3.txt"],
        outputs=["final_output.txt"],
        options={},
        working_dir="/some/dir",
    )
    graph = Graph.from_targets(
        {"target1": target1, "target2": target2, "target3": target3, "target4": target4}
    )
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    monkeypatch.setattr(scheduler, "should_run", lambda t: True)
    assert scheduler.schedule(target4)
    assert len(backend.submit.call_args_list) == 4
    assert call(target1, dependencies=set()) in backend.submit.call_args_list
    assert call(target2, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target3, dependencies=set([target2])) in backend.submit.call_args_list
    assert call(target4, dependencies=set([target3])) in backend.submit.call_args_list


def test_scheduling_branch_and_join_structure(backend, monkeypatch):
    target1 = Target(
        "TestTarget1",
        inputs=[],
        outputs=["output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    target2 = Target(
        "TestTarget2",
        inputs=["output1.txt"],
        outputs=["output2.txt"],
        options={},
        working_dir="/some/dir",
    )
    target3 = Target(
        "TestTarget3",
        inputs=["output1.txt"],
        outputs=["output3.txt"],
        options={},
        working_dir="/some/dir",
    )
    target4 = Target(
        "TestTarget4",
        inputs=["output2.txt", "output3.txt"],
        outputs=["final.txt"],
        options={},
        working_dir="/some/dir",
    )
    graph = Graph.from_targets(
        {"target1": target1, "target2": target2, "target3": target3, "target4": target4}
    )
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    monkeypatch.setattr(scheduler, "should_run", lambda t: True)
    assert scheduler.schedule(target4)
    assert len(backend.submit.call_args_list) == 4
    assert call(target1, dependencies=set([])) in backend.submit.call_args_list
    assert call(target2, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target3, dependencies=set([target1])) in backend.submit.call_args_list
    assert (
        call(target4, dependencies=set([target3, target2]))
        in backend.submit.call_args_list
    )


def test_scheduling_branch_and_join_structure_with_previously_submitted_dependency(
    backend, monkeypatch
):
    target1 = Target(
        "TestTarget1",
        inputs=[],
        outputs=["output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    target2 = Target(
        "TestTarget2",
        inputs=["output1.txt"],
        outputs=["output2.txt"],
        options={},
        working_dir="/some/dir",
    )
    target3 = Target(
        "TestTarget3",
        inputs=["output1.txt"],
        outputs=["output3.txt"],
        options={},
        working_dir="/some/dir",
    )
    target4 = Target(
        "TestTarget4",
        inputs=["output2.txt", "output3.txt"],
        outputs=["final.txt"],
        options={},
        working_dir="/some/dir",
    )

    graph = Graph.from_targets(
        {"target1": target1, "target2": target2, "target3": target3, "target4": target4}
    )
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    monkeypatch.setattr(scheduler, "should_run", lambda t: True)

    backend.submit(target1, dependencies=set())

    assert scheduler.schedule(target4)
    assert len(backend.submit.call_args_list) == 4
    assert call(target2, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target3, dependencies=set([target1])) in backend.submit.call_args_list
    assert (
        call(target4, dependencies=set([target3, target2]))
        in backend.submit.call_args_list
    )


def test_scheduling_non_submitted_targets_that_should_not_run(backend, monkeypatch):
    target1 = Target(
        "TestTarget1",
        inputs=[],
        outputs=["test_output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    target2 = Target(
        "TestTarget2",
        inputs=[],
        outputs=["test_output2.txt"],
        options={},
        working_dir="/some/dir",
    )
    target3 = Target(
        "TestTarget3",
        inputs=["test_output1.txt", "test_output2.txt"],
        outputs=["test_output3.txt"],
        options={},
        working_dir="/some/dir",
    )
    graph = Graph.from_targets(
        {"TestTarget1": target1, "TestTarget2": target2, "TestTarget3": target3}
    )
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    monkeypatch.setattr(scheduler, "should_run", lambda t: False)
    assert not scheduler.schedule(target3)
    assert backend.submit.call_args_list == []


def test_scheduling_many_targets_calls_schedule_for_each_target(backend, monkeypatch):
    target1 = Target(
        "TestTarget1",
        inputs=[],
        outputs=["test_output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    target2 = Target(
        "TestTarget2",
        inputs=[],
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
        inputs=["test_output2.txt"],
        outputs=["test_output4.txt"],
        options={},
        working_dir="/some/dir",
    )
    graph = Graph.from_targets(
        {
            "TestTarget1": target1,
            "TestTarget2": target2,
            "TestTarget3": target3,
            "TestTarget4": target4,
        }
    )
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, "should_run", lambda t: True)

    assert scheduler.schedule_many([target3, target4]) == [True, True]
    assert call(target4, dependencies=set([target2])) in backend.submit.call_args_list
    assert call(target3, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target2, dependencies=set()) in backend.submit.call_args_list
    assert call(target1, dependencies=set()) in backend.submit.call_args_list


def test_target_protected():
    target1 = Target(
        "TestTarget1",
        inputs=[],
        outputs=["test_output1.txt"],
        options={},
        working_dir="/some/dir",
        protect={"test_output1.txt"},
    )
    target2 = Target(
        "TestTarget2",
        inputs=[],
        outputs=["test_output2.txt"],
        options={},
        working_dir="/some/dir",
    )

    assert target1.protected == {"test_output1.txt"}
    assert target2.protected == set()


def test_map_naming_with_template_function():
    def my_template(path):
        return AnonymousTarget(inputs=[path], outputs=[path + ".new"], options={})

    files = ["a", "b", "c"]

    workflow = Workflow(working_dir="/some/dir")
    workflow.map(my_template, files)

    assert len(workflow.targets) == 3
    assert "my_template_0" in workflow.targets
    assert "my_template_1" in workflow.targets
    assert "my_template_2" in workflow.targets


def test_map_naming_with_template_class_instance():
    class MyTemplate:
        def __call__(self, path):
            return AnonymousTarget(inputs=[path], outputs=[path + ".new"], options={})

    files = ["a", "b", "c"]

    workflow = Workflow(working_dir="/some/dir")
    workflow.map(MyTemplate(), files)

    assert len(workflow.targets) == 3
    assert "MyTemplate_0" in workflow.targets
    assert "MyTemplate_1" in workflow.targets
    assert "MyTemplate_2" in workflow.targets


def test_map_naming_with_invalid_template_arg():
    files = ["a", "b", "c"]

    workflow = Workflow(working_dir="/some/dir")

    with pytest.raises(ValueError):
        workflow.map(42, files)


def test_map_with_custom_naming_function():
    def my_template(path):
        return AnonymousTarget(
            inputs={"path": path}, outputs={"path": path + ".new"}, options={}
        )

    files = ["a", "b", "c"]

    workflow = Workflow(working_dir="/some/dir")
    workflow.map(
        my_template, files, name=lambda i, t: "foo_{}".format(t.inputs["path"])
    )

    assert len(workflow.targets) == 3
    assert "foo_a" in workflow.targets
    assert "foo_b" in workflow.targets
    assert "foo_c" in workflow.targets


def test_map_with_custom_naming_string():
    def my_template(path):
        return AnonymousTarget(
            inputs={"path": path}, outputs={"path": path + ".new"}, options={}
        )

    files = ["a", "b", "c"]

    workflow = Workflow(working_dir="/some/dir")
    workflow.map(my_template, files, name="bar")

    assert len(workflow.targets) == 3
    assert "bar_0" in workflow.targets
    assert "bar_1" in workflow.targets
    assert "bar_2" in workflow.targets


@pytest.fixture
def mock_template(mocker,):
    mock_template = mocker.MagicMock()
    mock_template.__name__ = "mock_template"
    mock_template.return_value = AnonymousTarget(inputs=[], outputs=[], options={})
    return mock_template


def test_map_arg_passing_list_of_strings(mocker, mock_template):
    files = ["a", "b", "c"]

    workflow = Workflow(working_dir="/some/dir")
    workflow.map(mock_template, files)

    mock_template.assert_has_calls(
        [mocker.call("a"), mocker.call("b"), mocker.call("c")], any_order=True
    )


def test_map_arg_passing_list_of_tuples(mocker, mock_template):
    files = [("a", "/foo"), ("b", "/foo"), ("c", "/foo")]

    workflow = Workflow(working_dir="/some/dir")
    workflow.map(mock_template, files)

    mock_template.assert_has_calls(
        [mocker.call("a", "/foo"), mocker.call("b", "/foo"), mocker.call("c", "/foo")],
        any_order=True,
    )


def test_map_arg_passing_list_of_dicts(mocker, mock_template):
    files = [
        {"path": "a", "output_dir": "foo/"},
        {"path": "b", "output_dir": "foo/"},
        {"path": "c", "output_dir": "foo/"},
    ]

    workflow = Workflow(working_dir="/some/dir")
    workflow.map(mock_template, files)

    mock_template.assert_has_calls(
        [
            mocker.call(path="a", output_dir="foo/"),
            mocker.call(path="b", output_dir="foo/"),
            mocker.call(path="c", output_dir="foo/"),
        ],
        any_order=True,
    )


def test_map_arg_passing_list_of_dicts_with_extra(mocker, mock_template):
    files = [{"path": "a"}, {"path": "b"}, {"path": "c"}]

    workflow = Workflow(working_dir="/some/dir")
    workflow.map(mock_template, files, extra={"output_dir": "foo/"})

    mock_template.assert_has_calls(
        [
            mocker.call(path="a", output_dir="foo/"),
            mocker.call(path="b", output_dir="foo/"),
            mocker.call(path="c", output_dir="foo/"),
        ],
        any_order=True,
    )


def test_target_list():
    def my_template(path):
        return AnonymousTarget(
            inputs={"path": path}, outputs={"path": path + ".new"}, options={}
        )

    files = ["a", "b", "c"]

    workflow = Workflow(working_dir="/some/dir")
    target_list = workflow.map(my_template, files)

    assert len(target_list) == 3

    assert len(target_list.outputs) == 3
    assert target_list.outputs == [
        {"path": "a.new"},
        {"path": "b.new"},
        {"path": "c.new"},
    ]

    assert len(target_list.inputs) == 3
    assert target_list.inputs == [{"path": "a"}, {"path": "b"}, {"path": "c"}]


class FakeBackend(Backend):
    option_defaults = {
        "cores": 1,
        "memory": "1g",
    }

    def submit(self, target, dependencies):
        pass

    def status(self, target):
        return Status.UNKNOWN


def test_scheduler_injects_target_defaults_into_target_options_on_submit(mocker):
    target1 = Target(
        "TestTarget1", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    target2 = Target(
        "TestTarget2",
        inputs=[],
        outputs=[],
        options={"cores": 32},
        working_dir="/some/dir",
    )

    backend = FakeBackend()
    mocker.patch.object(backend, "submit", autospec=True)

    graph = Graph.from_targets({"TestTarget1": target1, "TestTarget2": target2})
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())

    scheduler.schedule(target1)
    assert target1.options == {"cores": 1, "memory": "1g"}

    scheduler.schedule(target2)
    assert target2.options == {"cores": 32, "memory": "1g"}


def test_scheduler_warns_user_when_submitting_target_with_unsupported_option(
    mocker, caplog
):
    target = Target(
        "TestTarget",
        inputs=[],
        outputs=[],
        options={"foo": "bar"},
        working_dir="/some/dir",
    )

    backend = FakeBackend()
    mocker.patch.object(backend, "submit", autospec=True)

    graph = Graph.from_targets({"TestTarget": target})

    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    scheduler.schedule(target)

    assert target.options == {"cores": 1, "memory": "1g"}
    assert caplog.record_tuples == [
        (
            "gwf.core",
            logging.WARNING,
            'Option "foo" used in "TestTarget" is not supported by backend. Ignored.',
        )
    ]


def test_scheduler_removes_options_with_none_value(mocker):
    target = Target(
        "TestTarget",
        inputs=[],
        outputs=[],
        options={"cores": None},
        working_dir="/some/dir",
    )

    backend = FakeBackend()
    mocker.patch.object(backend, "submit", autospec=True)

    graph = Graph.from_targets({"TestTarget": target})

    scheduler = Scheduler(graph=graph, backend=backend, filesystem=FakeFilesystem())
    scheduler.schedule(target)

    assert target.options == {"memory": "1g"}
