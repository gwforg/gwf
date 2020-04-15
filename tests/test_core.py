import logging
import unittest

import pytest

from gwf import AnonymousTarget, Graph, Target, Workflow
from gwf.exceptions import NameError, WorkflowError
from gwf.workflow import _flatten


@pytest.mark.parametrize(
    "value,expected",
    [
        (["a", "b"], ["a", "b"]),
        ({"A": ["a1", "a2"]}, ["a1", "a2"]),
        ({"A": ["a1", "a2"], "B": ["b1", "b2"]}, ["a1", "a2", "b1", "b2"]),
        ([{"A": "a1", "B": "b1"}, {"A": "a2", "B": "b2"}], ["a1", "a2", "b1", "b2"]),
        (
            [{"A": ["a1", "A1"], "B": "b1"}, {"A": ["a2", "A2"], "B": "b2"}],
            ["a1", "a2", "b1", "b2", "A1", "A2"],
        ),
    ],
)
def test_flatten(value, expected):
    assert set(_flatten(value)) == set(expected)


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


def test_graph_construction(graph_factory):
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

    graph = graph_factory([t1, t2, t3, t4])

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


def test_graph_raises_circular_dependency_error(graph_factory):
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


def test_graph_raises_when_two_targets_output_the_same_file(graph_factory):
    target1 = Target(
        "TestTarget1",
        inputs=[],
        outputs=["/some/dir/test_output.txt"],
        options={},
        working_dir="/some/dir",
    )
    target2 = Target(
        "TestTarget2",
        inputs=[],
        outputs=["test_output.txt"],
        options={},
        working_dir="/some/dir",
    )

    with pytest.raises(WorkflowError):
        graph_factory([target1, target2])


def test_graph_subset(diamond_graph):
    target1 = diamond_graph.targets["TestTarget1"]
    target2 = diamond_graph.targets["TestTarget2"]
    target3 = diamond_graph.targets["TestTarget3"]
    target4 = diamond_graph.targets["TestTarget4"]

    g1 = diamond_graph.subset([target3])
    assert set([target1, target3]) == set(g1.targets.values())

    g2 = diamond_graph.subset([target2])
    assert set([target1, target2]) == set(g2.targets.values())

    g3 = diamond_graph.subset([target4])
    assert set([target1, target2, target3, target4]) == set(g3.targets.values())


def test_schedule_if_one_of_its_output_files_does_not_exist(diamond_graph, schedule):
    target = diamond_graph.targets["TestTarget1"]
    scheduled, reasons = schedule([target], diamond_graph)
    assert target in scheduled
    assert scheduled[target] == set()
    assert reasons[target] == (
        "TestTarget1 was scheduled because its output "
        "file /some/dir/test_output1.txt does not exist"
    )


def test_schedule_if_one_of_its_dependencies_was_scheduled(diamond_graph, schedule):
    target1 = diamond_graph.targets["TestTarget1"]
    target2 = diamond_graph.targets["TestTarget2"]
    scheduled, reasons = schedule([target2], diamond_graph)

    assert target2 in scheduled
    assert reasons[target2] == (
        "TestTarget2 was scheduled because its dependency TestTarget1 was scheduled"
    )
    assert scheduled[target2] == set([target1])


def test_schedule_if_it_is_a_sink(trivial_graph, schedule):
    target = trivial_graph.targets["TestTarget"]
    scheduled, reasons = schedule([target], trivial_graph)

    assert target in scheduled
    assert reasons[target] == "TestTarget was scheduled because it is a sink"
    assert scheduled[target] == set()


def test_schedule_if_any_input_file_is_newer_than_any_output_file(
    diamond_graph, filesystem, schedule
):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=0)
    filesystem.add_file("/some/dir/test_output2.txt", changed_at=1)
    filesystem.add_file("/some/dir/test_output3.txt", changed_at=3)
    filesystem.add_file("/some/dir/final_output.txt", changed_at=2)

    target4 = diamond_graph.targets["TestTarget4"]

    scheduled, reasons = schedule([target4], diamond_graph)
    assert len(scheduled) == 1

    assert target4 in scheduled
    assert scheduled[target4] == set()
    assert reasons[target4] == (
        "TestTarget4 was scheduled because input file /some/dir/test_output3.txt "
        "is newer than output file /some/dir/final_output.txt"
    )


def test_do_not_schedule_if_it_is_a_source_and_all_outputs_exist(
    diamond_graph, filesystem, schedule
):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=1)

    target = diamond_graph.targets["TestTarget1"]
    scheduled, reasons = schedule([target], diamond_graph)
    assert target not in scheduled
    assert reasons[target] == "TestTarget1 was not scheduled because it is a source"


def test_do_not_schedule_if_all_outputs_are_newer_then_the_inputs(
    diamond_graph, schedule, filesystem
):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=0)
    filesystem.add_file("/some/dir/test_output2.txt", changed_at=1)
    filesystem.add_file("/some/dir/test_output3.txt", changed_at=3)
    filesystem.add_file("/some/dir/final_output.txt", changed_at=4)

    scheduled, not_scheduled = schedule(
        [diamond_graph.targets["TestTarget4"]], diamond_graph
    )
    assert len(scheduled) == 0
    assert len(not_scheduled) == 4


def test_scheduler_raises_if_input_file_is_not_provided_and_does_not_exist(
    schedule, graph_factory
):
    target = Target(
        "TestTarget",
        inputs=["in.txt"],
        outputs=["out.txt"],
        options={},
        working_dir="/some/dir",
    )
    graph = graph_factory([target])

    with pytest.raises(WorkflowError):
        schedule([target], graph)


def test_scheduling_target_with_deep_deps_that_are_not_submitted(
    schedule, graph_factory
):
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

    graph = graph_factory([target1, target2, target3, target4])
    scheduled, reasons = schedule([target4], graph)

    assert len(scheduled) == 4

    assert scheduled[target1] == set([])
    assert scheduled[target2] == set([target1])
    assert scheduled[target3] == set([target2])
    assert scheduled[target4] == set([target3])


def test_scheduling_branch_and_join_structure(schedule, graph_factory):
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
    graph = graph_factory([target1, target2, target3, target4])
    scheduled, reasons = schedule([target4], graph)

    assert len(scheduled) == 4

    assert scheduled[target1] == set()
    assert scheduled[target2] == set([target1])
    assert scheduled[target3] == set([target1])
    assert scheduled[target4] == set([target2, target3])


def test_scheduling_multiple_targets(diamond_graph, schedule):
    scheduled, reasons = schedule(
        [diamond_graph.targets["TestTarget3"], diamond_graph.targets["TestTarget2"]],
        diamond_graph,
    )
    assert len(scheduled) == 3


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


@pytest.mark.skip(msg="injection of backend defaults will be moved to backends")
def test_scheduler_injects_target_defaults_into_target_options_on_submit(
    backend, filesystem
):
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

    graph = Graph.from_targets({"TestTarget1": target1, "TestTarget2": target2})
    scheduler = Scheduler(graph=graph, backend=backend, filesystem=filesystem)

    scheduler.schedule(target1)
    assert target1.options == {"cores": 1, "memory": "1g"}

    scheduler.schedule(target2)
    assert target2.options == {"cores": 32, "memory": "1g"}


@pytest.mark.skip(msg="injection of backend defaults will be moved to backends")
def test_scheduler_warns_user_when_submitting_target_with_unsupported_option(
    backend, caplog, filesystem
):
    target = Target(
        "TestTarget",
        inputs=[],
        outputs=[],
        options={"foo": "bar"},
        working_dir="/some/dir",
    )

    graph = Graph.from_targets({"TestTarget": target})

    scheduler = Scheduler(graph=graph, backend=backend, filesystem=filesystem)
    scheduler.schedule(target)

    assert target.options == {"cores": 1, "memory": "1g"}
    assert caplog.record_tuples == [
        (
            "gwf.core",
            logging.WARNING,
            'Option "foo" used in "TestTarget" is not supported by backend. Ignored.',
        )
    ]


@pytest.mark.skip(msg="injection of backend defaults will be moved to backends")
def test_scheduler_removes_options_with_none_value(backend, filesystem):
    target = Target(
        "TestTarget",
        inputs=[],
        outputs=[],
        options={"cores": None},
        working_dir="/some/dir",
    )

    graph = Graph.from_targets({"TestTarget": target})

    scheduler = Scheduler(graph=graph, backend=backend, filesystem=filesystem)
    scheduler.schedule(target)

    assert target.options == {"memory": "1g"}
