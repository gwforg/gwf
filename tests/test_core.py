import unittest

import pytest

from gwf.core import (
    CircularDependencyError,
    FileProvidedByMultipleTargetsError,
    Graph,
    Target,
    _flatten,
)
from gwf.exceptions import NameError, WorkflowError
from gwf.scheduling import Reason


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

    def test_target_with_protected_outputs(self):
        target = Target(
            name="TestTarget",
            inputs=["test_input1.txt", "test_input2.txt"],
            outputs=[],
            options={},
            working_dir="/some/path",
            protect=["test_input1.txt"],
        )
        self.assertEqual(target.protected(), set(["/some/path/test_input1.txt"]))

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

    with pytest.raises(FileProvidedByMultipleTargetsError):
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
    with pytest.raises(CircularDependencyError):
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

    with pytest.raises(FileProvidedByMultipleTargetsError):
        graph_factory([target1, target2])


def test_schedule_if_one_of_its_output_files_does_not_exist(
    diamond_graph, filesystem, schedule
):
    target = diamond_graph.targets["TestTarget1"]
    reasons = schedule(diamond_graph, endpoints=[target], fs=filesystem)
    reason = reasons[target]
    assert isinstance(reason, Reason.MissingOutput)
    assert reason.scheduled


def test_schedule_if_one_of_its_dependencies_was_scheduled(diamond_graph, schedule):
    target2 = diamond_graph.targets["TestTarget2"]
    reasons = schedule(diamond_graph, endpoints=[target2])
    reason = reasons[target2]
    assert isinstance(reason, Reason.DependencyScheduled)
    assert reason.target is target2
    assert reason.scheduled


@pytest.mark.skip()
def test_schedule_if_it_is_a_sink(trivial_graph, schedule):
    target = trivial_graph.targets["TestTarget"]
    reasons = schedule(trivial_graph, endpoints=[target])
    reason = reasons[target]
    assert reason.scheduled
    assert isinstance(reason, Reason.IsSink)


def test_schedule_if_any_input_file_is_newer_than_any_output_file(
    diamond_graph, filesystem, schedule
):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=0)
    filesystem.add_file("/some/dir/test_output2.txt", changed_at=1)
    filesystem.add_file("/some/dir/test_output3.txt", changed_at=3)
    filesystem.add_file("/some/dir/final_output.txt", changed_at=2)

    target4 = diamond_graph.targets["TestTarget4"]
    reasons = schedule(diamond_graph, endpoints=[target4])
    reason = reasons[target4]
    assert reason.target is target4
    assert reason.scheduled
    assert isinstance(reason, Reason.OutOfDate)


def test_do_not_schedule_if_it_is_a_source_and_all_outputs_exist(
    diamond_graph, filesystem, schedule
):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=1)

    target = diamond_graph.targets["TestTarget1"]
    reasons = schedule(diamond_graph, endpoints=[target])
    reason = reasons[target]
    assert not reason.scheduled
    assert isinstance(reason, Reason.IsSource)


def test_do_not_schedule_if_all_outputs_are_newer_then_the_inputs(
    diamond_graph, schedule, filesystem
):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=0)
    filesystem.add_file("/some/dir/test_output2.txt", changed_at=1)
    filesystem.add_file("/some/dir/test_output3.txt", changed_at=3)
    filesystem.add_file("/some/dir/final_output.txt", changed_at=4)

    target = diamond_graph.targets["TestTarget4"]
    reasons = schedule(diamond_graph, endpoints=[target])
    reason = reasons[target]
    assert not reason.scheduled


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
        schedule(graph, endpoints=[target])


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
    reasons = schedule(graph, endpoints=[target4])

    assert reasons[target1].scheduled
    assert reasons[target2].scheduled
    assert reasons[target3].scheduled
    assert reasons[target4].scheduled


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
    reasons = schedule(graph, endpoints=[target4])

    assert reasons[target1].scheduled
    assert reasons[target2].scheduled
    assert reasons[target3].scheduled
    assert reasons[target4].scheduled


@pytest.mark.skip()
def test_scheduling_multiple_targets(diamond_graph, schedule):
    reasons = schedule(
        diamond_graph,
        endpoints=[
            diamond_graph.targets["TestTarget3"],
            diamond_graph.targets["TestTarget2"],
        ],
    )
    # assert len(reasons) == 2


def test_scheduling_when_spec_changed(diamond_graph, schedule, filesystem):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=1)
    target = diamond_graph.targets["TestTarget1"]
    reasons = schedule(diamond_graph, endpoints=[target], spec_hashes={})
    reason = reasons[target]
    assert reason.scheduled
    assert isinstance(reason, Reason.SpecChanged)


def test_scheduling_when_spec_not_changed(diamond_graph, schedule, filesystem):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=1)
    target = diamond_graph.targets["TestTarget1"]
    reasons = schedule(
        diamond_graph,
        endpoints=[target],
        spec_hashes={target.name: "da39a3ee5e6b4b0d3255bfef95601890afd80709"},
    )
    reason = reasons[target]
    assert not reason.scheduled
