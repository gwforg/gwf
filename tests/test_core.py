import unittest

import pytest

from gwf.core import (
    CircularDependencyError,
    FileProvidedByMultipleTargetsError,
    Graph,
    InvalidPathError,
    Module,
    NoopSpecHashes,
    Status,
    Target,
    TempFileUsedOutsideModuleError,
    UnresolvedInputError,
    _flatten,
)
from gwf.exceptions import GWFError
from gwf.scheduling import get_status_map
from gwf.temp import _clear_temp_registry, is_temp, temp


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
        with self.assertRaises(GWFError):
            Target(
                "123abc", inputs=[], outputs=[], options={}, working_dir="/some/path"
            )

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

    def test_str_on_target(self):
        target = Target(
            "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/path"
        )
        self.assertEqual(str(target), "TestTarget")

    def test_input_is_empty_string(self):
        with self.assertRaises(InvalidPathError):
            Target(
                name="TestTarget",
                inputs=["ab", ""],
                outputs=[],
                options={},
                working_dir="/some/path",
            )

    def test_output_is_empty_string(self):
        with self.assertRaises(InvalidPathError):
            Target(
                name="TestTarget",
                inputs=[],
                outputs=["ab", ""],
                options={},
                working_dir="/some/path",
            )

    def test_input_contains_nonprintable(self):
        with self.assertRaises(InvalidPathError):
            Target(
                name="TestTarget",
                inputs=["a\nb", "ac"],
                outputs=[],
                options={},
                working_dir="/some/path",
            )

    def test_output_contains_nonprintable(self):
        with self.assertRaises(InvalidPathError):
            Target(
                name="TestTarget",
                inputs=[],
                outputs=["a\nb", "ac"],
                options={},
                working_dir="/some/path",
            )


class TestModule(unittest.TestCase):
    def test_module_with_invalid_name_raises_exception(self):
        with self.assertRaises(GWFError):
            Module("123abc", targets=[])

    def test_module_flattened_outputs(self):
        target1 = Target(
            name="Target1",
            inputs=[],
            outputs=["out1.txt"],
            options={},
            working_dir="/some/dir",
        )
        target2 = Target(
            name="Target2",
            inputs=[],
            outputs=["out2.txt"],
            options={},
            working_dir="/some/dir",
        )
        module = Module(name="TestModule", targets=[target1, target2])

        outputs = module.flattened_outputs()
        assert "/some/dir/out1.txt" in outputs
        assert "/some/dir/out2.txt" in outputs

    def test_module_flattened_outputs_nested(self):
        target1 = Target(
            name="Target1",
            inputs=[],
            outputs=["out1.txt"],
            options={},
            working_dir="/some/dir",
        )
        target2 = Target(
            name="Target2",
            inputs=[],
            outputs=["out2.txt"],
            options={},
            working_dir="/some/dir",
        )
        inner_module = Module(name="InnerModule", targets=[target2])
        outer_module = Module(name="OuterModule", targets=[target1, inner_module])

        outputs = outer_module.flattened_outputs()
        assert "/some/dir/out1.txt" in outputs
        assert "/some/dir/out2.txt" in outputs

    def test_module_all_targets(self):
        target1 = Target(
            name="Target1",
            inputs=[],
            outputs=["out1.txt"],
            options={},
            working_dir="/some/dir",
        )
        target2 = Target(
            name="Target2",
            inputs=[],
            outputs=["out2.txt"],
            options={},
            working_dir="/some/dir",
        )
        inner_module = Module(name="InnerModule", targets=[target2])
        outer_module = Module(name="OuterModule", targets=[target1, inner_module])

        all_targets = outer_module.all_targets()
        assert target1 in all_targets
        assert target2 in all_targets
        assert len(all_targets) == 2


class TestTempFunction:
    def test_temp_marks_path_as_temporary(self):
        path = "some/path.txt"
        result = temp(path)

        assert result is path
        assert is_temp(result) is True

        _clear_temp_registry()

    def test_regular_path_is_not_temp(self):
        path = "some/path.txt"

        assert is_temp(path) is False

        _clear_temp_registry()


class TestTempFileBoundaryValidation:
    def test_temp_chain_within_module_allowed(self, filesystem):
        step1 = Target(
            name="Step1",
            inputs=[],
            outputs=[temp("temp1.txt")],
            options={},
            working_dir="/some/dir",
        )
        step2 = Target(
            name="Step2",
            inputs=["temp1.txt"],
            outputs=[temp("temp2.txt")],
            options={},
            working_dir="/some/dir",
        )
        step3 = Target(
            name="Step3",
            inputs=["temp2.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )
        module = Module(name="Pipeline", targets=[step1, step2, step3])

        graph = Graph.from_targets([module], filesystem)
        assert len(graph.targets) == 3

        _clear_temp_registry()

    def test_multiple_consumers_same_temp_within_module_allowed(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("shared_temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        consumer1 = Target(
            name="Consumer1",
            inputs=["shared_temp.txt"],
            outputs=["output1.txt"],
            options={},
            working_dir="/some/dir",
        )
        consumer2 = Target(
            name="Consumer2",
            inputs=["shared_temp.txt"],
            outputs=["output2.txt"],
            options={},
            working_dir="/some/dir",
        )
        module = Module(name="TestModule", targets=[producer, consumer1, consumer2])

        graph = Graph.from_targets([module], filesystem)
        assert len(graph.targets) == 3

        _clear_temp_registry()

    def test_multiple_consumers_one_outside_module_raises(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("shared_temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        consumer_inside = Target(
            name="ConsumerInside",
            inputs=["shared_temp.txt"],
            outputs=["output1.txt"],
            options={},
            working_dir="/some/dir",
        )
        module = Module(name="TestModule", targets=[producer, consumer_inside])

        consumer_outside = Target(
            name="ConsumerOutside",
            inputs=["shared_temp.txt"],
            outputs=["output2.txt"],
            options={},
            working_dir="/some/dir",
        )

        with pytest.raises(TempFileUsedOutsideModuleError):
            Graph.from_targets([module, consumer_outside], filesystem)

        _clear_temp_registry()

    def test_diamond_dependency_with_temp_in_module_allowed(self, filesystem):
        root = Target(
            name="Root",
            inputs=[],
            outputs=[temp("root_out.txt")],
            options={},
            working_dir="/some/dir",
        )
        left = Target(
            name="Left",
            inputs=["root_out.txt"],
            outputs=[temp("left_out.txt")],
            options={},
            working_dir="/some/dir",
        )
        right = Target(
            name="Right",
            inputs=["root_out.txt"],
            outputs=[temp("right_out.txt")],
            options={},
            working_dir="/some/dir",
        )
        join = Target(
            name="Join",
            inputs=["left_out.txt", "right_out.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )
        module = Module(name="Diamond", targets=[root, left, right, join])

        graph = Graph.from_targets([module], filesystem)
        assert len(graph.targets) == 4

        _clear_temp_registry()

    def test_producer_in_inner_module_consumer_in_outer_allowed(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("intermediate.txt")],
            options={},
            working_dir="/some/dir",
        )
        inner_module = Module(name="InnerModule", targets=[producer])

        consumer = Target(
            name="Consumer",
            inputs=["intermediate.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )
        outer_module = Module(name="OuterModule", targets=[inner_module, consumer])

        graph = Graph.from_targets([outer_module], filesystem)
        assert len(graph.targets) == 2

        _clear_temp_registry()

    def test_deeply_nested_modules_temp_allowed(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("deep_temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        consumer = Target(
            name="Consumer",
            inputs=["deep_temp.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )

        level3 = Module(name="Level3", targets=[consumer])
        level2 = Module(name="Level2", targets=[level3])
        level1 = Module(name="Level1", targets=[producer, level2])

        graph = Graph.from_targets([level1], filesystem)
        assert len(graph.targets) == 2

        _clear_temp_registry()

    def test_separate_branches_under_parent_allowed(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("branch_temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        branch_a = Module(name="BranchA", targets=[producer])

        consumer = Target(
            name="Consumer",
            inputs=["branch_temp.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )
        branch_b = Module(name="BranchB", targets=[consumer])

        parent = Module(name="Parent", targets=[branch_a, branch_b])

        graph = Graph.from_targets([parent], filesystem)
        assert len(graph.targets) == 2

        _clear_temp_registry()

    def test_multiple_temp_files_mixed_validity(self, filesystem):
        producer1 = Target(
            name="Producer1",
            inputs=[],
            outputs=[temp("temp1.txt")],
            options={},
            working_dir="/some/dir",
        )
        consumer1 = Target(
            name="Consumer1",
            inputs=["temp1.txt"],
            outputs=["out1.txt"],
            options={},
            working_dir="/some/dir",
        )
        module_a = Module(name="ModuleA", targets=[producer1, consumer1])

        producer2 = Target(
            name="Producer2",
            inputs=[],
            outputs=[temp("temp2.txt")],
            options={},
            working_dir="/some/dir",
        )
        module_b = Module(name="ModuleB", targets=[producer2])

        consumer2 = Target(
            name="Consumer2",
            inputs=["temp2.txt"],
            outputs=["out2.txt"],
            options={},
            working_dir="/some/dir",
        )

        with pytest.raises(TempFileUsedOutsideModuleError):
            Graph.from_targets([module_a, module_b, consumer2], filesystem)

        _clear_temp_registry()

    def test_global_temp_global_consumer_not_allowed(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("global_temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        consumer = Target(
            name="Consumer",
            inputs=["global_temp.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )

        with pytest.raises(TempFileUsedOutsideModuleError):
            Graph.from_targets([producer, consumer], filesystem)

        _clear_temp_registry()

    def test_temp_with_different_working_dirs(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("output.txt")],
            options={},
            working_dir="/dir_a",
        )
        module = Module(name="TestModule", targets=[producer])

        consumer = Target(
            name="Consumer",
            inputs=["output.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/dir_b",
        )
        consumer_module = Module(name="ConsumerModule", targets=[consumer])

        filesystem.add_file("/dir_b/output.txt", changed_at=1)

        graph = Graph.from_targets([module, consumer_module], filesystem)
        assert len(graph.targets) == 2

        _clear_temp_registry()

    def test_empty_module_no_error(self, filesystem):
        empty_module = Module(name="EmptyModule", targets=[])
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        consumer = Target(
            name="Consumer",
            inputs=["temp.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )
        # Temp files must be in a module
        work_module = Module(name="WorkModule", targets=[producer, consumer])

        graph = Graph.from_targets([empty_module, work_module], filesystem)
        assert len(graph.targets) == 2

        _clear_temp_registry()

    def test_only_regular_files_no_validation_needed(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=["regular.txt"],
            options={},
            working_dir="/some/dir",
        )
        module_a = Module(name="ModuleA", targets=[producer])

        consumer = Target(
            name="Consumer",
            inputs=["regular.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )
        module_b = Module(name="ModuleB", targets=[consumer])

        graph = Graph.from_targets([module_a, module_b], filesystem)
        assert len(graph.targets) == 2

        _clear_temp_registry()

    def test_consumer_module_has_own_producer_allowed(self, filesystem):
        global_producer = Target(
            name="GlobalProducer",
            inputs=[],
            outputs=[temp("shared_temp.txt")],
            options={},
            working_dir="/some/dir",
        )

        module_producer = Target(
            name="GlobalProducer",
            inputs=[],
            outputs=[temp("shared_temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        consumer = Target(
            name="Consumer",
            inputs=["shared_temp.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )
        module = Module(name="TestModule", targets=[module_producer, consumer])

        graph = Graph.from_targets([global_producer, module], filesystem)
        assert len(graph.targets) == 2

        _clear_temp_registry()

    def test_temp_file_module_to_global_raises(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("module_temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        module = Module(name="TestModule", targets=[producer])

        consumer = Target(
            name="Consumer",
            inputs=["module_temp.txt"],
            outputs=["final.txt"],
            options={},
            working_dir="/some/dir",
        )

        with pytest.raises(TempFileUsedOutsideModuleError):
            Graph.from_targets([module, consumer], filesystem)

        _clear_temp_registry()

    def test_global_producer_temp_not_allowed(self, filesystem):
        producer = Target(
            name="Producer",
            inputs=[],
            outputs=[temp("global_temp.txt")],
            options={},
            working_dir="/some/dir",
        )

        with pytest.raises(TempFileUsedOutsideModuleError):
            Graph.from_targets([producer], filesystem)

        _clear_temp_registry()

    def test_complex_workflow_with_multiple_modules_and_temps(self, filesystem):
        a_step1 = Target(
            name="A_Step1",
            inputs=[],
            outputs=[temp("a_temp1.txt"), "a_output1.txt"],
            options={},
            working_dir="/some/dir",
        )
        a_step2 = Target(
            name="A_Step2",
            inputs=["a_temp1.txt"],
            outputs=["a_final.txt"],
            options={},
            working_dir="/some/dir",
        )
        module_a = Module(name="ModuleA", targets=[a_step1, a_step2])

        b_step1 = Target(
            name="B_Step1",
            inputs=["a_output1.txt"],
            outputs=[temp("b_temp.txt")],
            options={},
            working_dir="/some/dir",
        )
        b_step2 = Target(
            name="B_Step2",
            inputs=["b_temp.txt"],
            outputs=["b_final.txt"],
            options={},
            working_dir="/some/dir",
        )
        module_b = Module(name="ModuleB", targets=[b_step1, b_step2])

        final = Target(
            name="Final",
            inputs=["a_final.txt", "b_final.txt"],
            outputs=["workflow_complete.txt"],
            options={},
            working_dir="/some/dir",
        )

        graph = Graph.from_targets([module_a, module_b, final], filesystem)
        assert len(graph.targets) == 5

        _clear_temp_registry()


def test_module_is_complete(filesystem):
    target1 = Target(
        name="Target1",
        inputs=[],
        outputs=["out1.txt"],
        options={},
        working_dir="/some/dir",
    )
    module = Module(name="TestModule", targets=[target1])

    assert module.is_complete(filesystem) is False

    filesystem.add_file("/some/dir/out1.txt", changed_at=1)
    assert module.is_complete(filesystem) is True


def test_graph_construction(filesystem):
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

    filesystem.add_file("/some/dir/test_input1.txt", 0)
    filesystem.add_file("/some/dir/test_input2.txt", 0)

    graph = Graph.from_targets([t1, t2, t3, t4], filesystem)

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


def test_graph_construction_with_modules(filesystem):
    target1 = Target(
        name="Target1",
        inputs=["t1_input1.txt"],
        outputs=["t1_output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    target2 = Target(
        name="Target2",
        inputs=["t1_output1.txt"],
        outputs=["t2_output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    target3 = Target(
        name="Target3",
        inputs=["t2_output1.txt"],
        outputs=["t3_output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    module2 = Module(
        name="Module2",
        targets=[target3],
    )
    module1 = Module(
        name="Module1",
        targets=[target2, module2],
    )

    filesystem.add_file("/some/dir/t1_input1.txt", changed_at=0)
    graph = Graph.from_targets([target1, module1], filesystem)

    assert len(graph.targets) == 3
    assert target1 in graph.targets.values()
    assert target2 in graph.targets.values()
    assert target3 in graph.targets.values()

    assert graph.dependencies[target1] == set()
    assert graph.dependents[target1] == {target2}

    assert graph.dependencies[target2] == {target1}
    assert graph.dependents[target2] == {target3}

    assert graph.dependencies[target3] == {target2}
    assert graph.dependents[target3] == set()

    assert "Module1" in graph.modules
    assert "Module2" in graph.modules

    assert "Module1" in graph.target_modules["Target2"]
    assert "Module1" in graph.target_modules["Target3"]
    assert "Module2" in graph.target_modules["Target3"]


def test_graph_construction_with_completed_module(filesystem):
    target1 = Target(
        name="Target1",
        inputs=["t1_input1.txt"],
        outputs=["t1_output1.txt"],
        options={},
        working_dir="/some/dir",
    )
    module1 = Module(
        name="Module1",
        targets=[target1],
    )
    target2 = Target(
        name="Target2",
        inputs=["t1_output1.txt"],
        outputs=["t2_output1.txt"],
        options={},
        working_dir="/some/dir",
    )

    filesystem.add_file("/some/dir/t1_input1.txt", changed_at=0)
    filesystem.add_file("/some/dir/t1_output1.txt", changed_at=1)
    graph = Graph.from_targets([module1, target2], filesystem)

    assert len(graph.targets) == 2
    assert target1 in graph.targets.values()
    assert target2 in graph.targets.values()

    assert graph.dependencies[target1] == set()
    assert graph.dependents[target1] == {target2}

    assert graph.dependencies[target2] == {target1}
    assert graph.dependents[target2] == set()

    assert "Module1" in graph.modules
    assert graph.modules["Module1"].is_complete(filesystem)

    assert "Module1" in graph.target_modules["Target1"]


def test_graph_target_in_multiple_modules(filesystem):
    shared_target = Target(
        name="SharedTarget",
        inputs=["input.txt"],
        outputs=["shared_output.txt"],
        options={},
        working_dir="/some/dir",
    )
    target_a = Target(
        name="TargetA",
        inputs=["shared_output.txt"],
        outputs=["output_a.txt"],
        options={},
        working_dir="/some/dir",
    )
    target_b = Target(
        name="TargetB",
        inputs=["shared_output.txt"],
        outputs=["output_b.txt"],
        options={},
        working_dir="/some/dir",
    )

    module_a = Module(name="ModuleA", targets=[shared_target, target_a])
    module_b = Module(name="ModuleB", targets=[shared_target, target_b])

    filesystem.add_file("/some/dir/input.txt", changed_at=0)
    graph = Graph.from_targets([module_a, module_b], filesystem)

    assert len(graph.targets) == 3
    assert shared_target in graph.targets.values()
    assert target_a in graph.targets.values()
    assert target_b in graph.targets.values()

    assert "ModuleA" in graph.target_modules["SharedTarget"]
    assert "ModuleB" in graph.target_modules["SharedTarget"]

    assert graph.target_modules["TargetA"] == {"ModuleA"}
    assert graph.target_modules["TargetB"] == {"ModuleB"}


def test_graph_get_modules(filesystem):
    target1 = Target(
        name="Target1",
        inputs=[],
        outputs=["out1.txt"],
        options={},
        working_dir="/some/dir",
    )
    module1 = Module(name="Module1", targets=[target1])
    module2 = Module(name="Module2", targets=[target1])

    graph = Graph.from_targets([module1, module2], filesystem)

    modules = graph.get_modules(target1)
    assert len(modules) == 2
    module_names = {m.name for m in modules}
    assert "Module1" in module_names
    assert "Module2" in module_names


def test_graph_raises_multiple_providers_error(filesystem):
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
        Graph.from_targets({"Target1": t1, "Target2": t2}, filesystem)


def test_graph_raises_circular_dependency_error(filesystem):
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
        Graph.from_targets({"Target1": t1, "Target2": t2, "Target3": t3}, filesystem)


def test_graph_raises_when_two_targets_output_the_same_file(filesystem):
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
        Graph.from_targets([target1, target2], filesystem)


def test_schedule_if_one_of_its_output_files_does_not_exist(
    diamond_graph,
    filesystem,
    backend,
):
    target = diamond_graph.targets["TestTarget1"]
    target_states = get_status_map(
        graph=diamond_graph,
        endpoints=[target],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target] == Status.SHOULDRUN


def test_schedule_if_one_of_its_dependencies_was_scheduled(
    diamond_graph, filesystem, backend
):
    target = diamond_graph.targets["TestTarget2"]
    target_states = get_status_map(
        graph=diamond_graph,
        endpoints=[target],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target] == Status.SHOULDRUN


def test_schedule_if_it_is_a_sink(trivial_graph, filesystem, backend):
    target = trivial_graph.targets["TestTarget"]
    target_states = get_status_map(
        graph=trivial_graph,
        endpoints=[target],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target] == Status.SHOULDRUN


def test_schedule_if_any_input_file_is_newer_than_any_output_file(
    diamond_graph, filesystem, backend
):
    filesystem.add_file("/some/dir/test_output1.txt", changed_at=0)
    filesystem.add_file("/some/dir/test_output2.txt", changed_at=1)
    filesystem.add_file("/some/dir/test_output3.txt", changed_at=3)
    filesystem.add_file("/some/dir/final_output.txt", changed_at=2)

    target = diamond_graph.targets["TestTarget4"]
    target_states = get_status_map(
        graph=diamond_graph,
        endpoints=[target],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target] == Status.SHOULDRUN


def test_schedule_if_it_is_a_source_and_has_missing_output_file(filesystem, backend):
    target = Target(
        name="Foo", inputs=[], outputs=["foo"], options={}, working_dir="/some/dir"
    )
    graph = Graph.from_targets([target], filesystem)

    target_states = get_status_map(
        graph=graph,
        endpoints=[target],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target] == Status.SHOULDRUN

    filesystem.add_file("/some/dir/foo", changed_at=1)
    target_states = get_status_map(
        graph=graph,
        endpoints=[target],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target] == Status.COMPLETED


def test_do_not_schedule_if_all_outputs_are_newer_then_the_inputs(backend, filesystem):
    filesystem.add_file("/some/dir/input1", changed_at=0)
    filesystem.add_file("/some/dir/input2", changed_at=1)
    filesystem.add_file("/some/dir/output1", changed_at=3)
    filesystem.add_file("/some/dir/output2", changed_at=4)

    target = Target(
        name="Foo",
        inputs=["input1", "input2"],
        outputs=["output1", "output2"],
        options={},
        working_dir="/some/dir",
    )
    graph = Graph.from_targets([target], filesystem)

    target_states = get_status_map(
        graph=graph,
        endpoints=[target],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target] == Status.COMPLETED


def test_graph_raises_if_input_file_is_not_provided_and_does_not_exist(filesystem):
    target = Target(
        "TestTarget",
        inputs=["in.txt"],
        outputs=["out.txt"],
        options={},
        working_dir="/some/dir",
    )

    with pytest.raises(UnresolvedInputError):
        Graph.from_targets([target], filesystem)


def test_scheduling_target_with_deep_deps_that_are_not_submitted(filesystem, backend):
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

    graph = Graph.from_targets([target1, target2, target3, target4], filesystem)
    target_states = get_status_map(
        graph=graph,
        endpoints=[target4],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target1] == Status.SHOULDRUN
    assert target_states[target2] == Status.SHOULDRUN
    assert target_states[target3] == Status.SHOULDRUN
    assert target_states[target4] == Status.SHOULDRUN


def test_scheduling_branch_and_join_structure(filesystem, backend):
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
    graph = Graph.from_targets([target1, target2, target3, target4], filesystem)
    target_states = get_status_map(
        graph=graph,
        endpoints=[target4],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )
    assert target_states[target1] == Status.SHOULDRUN
    assert target_states[target2] == Status.SHOULDRUN
    assert target_states[target3] == Status.SHOULDRUN
    assert target_states[target4] == Status.SHOULDRUN


def test_scheduling_with_spec_hashing(backend, spec_hashes, filesystem):
    filesystem.add_file("/some/dir/input.txt", changed_at=1)
    filesystem.add_file("/some/dir/output.txt", changed_at=2)

    target = Target(
        "TestTarget",
        inputs=["input.txt"],
        outputs=["output.txt"],
        options={},
        working_dir="/some/dir",
        spec="foo",
    )
    graph = Graph.from_targets([target], filesystem)

    target_states = get_status_map(
        graph=graph,
        fs=filesystem,
        backend=backend,
        spec_hashes=spec_hashes,
    )
    assert target_states[target] == Status.SHOULDRUN

    spec_hashes.update(target)

    target_states = get_status_map(
        graph=graph,
        fs=filesystem,
        backend=backend,
        spec_hashes=spec_hashes,
    )
    assert target_states[target] == Status.COMPLETED

    target.spec = "bar"

    target_states = get_status_map(
        graph=graph,
        fs=filesystem,
        backend=backend,
        spec_hashes=spec_hashes,
    )
    assert target_states[target] == Status.SHOULDRUN


def test_scheduling_skips_target_in_complete_module(filesystem, backend):
    """Target in a complete module should not be scheduled."""
    target1 = Target(
        name="Target1",
        inputs=[],
        outputs=["out1.txt"],
        options={},
        working_dir="/some/dir",
    )
    module1 = Module(name="Module1", targets=[target1])

    filesystem.add_file("/some/dir/out1.txt", changed_at=1)

    graph = Graph.from_targets([module1], filesystem)
    target_states = get_status_map(
        graph=graph,
        endpoints=[target1],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )

    assert target_states[target1] == Status.COMPLETED


def test_scheduling_runs_target_in_incomplete_module(filesystem, backend):
    """Target in an incomplete module should be scheduled."""
    target1 = Target(
        name="Target1",
        inputs=[],
        outputs=["out1.txt"],
        options={},
        working_dir="/some/dir",
    )
    module1 = Module(name="Module1", targets=[target1])

    graph = Graph.from_targets([module1], filesystem)
    target_states = get_status_map(
        graph=graph,
        endpoints=[target1],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )

    assert target_states[target1] == Status.SHOULDRUN


def test_scheduling_target_in_multiple_modules_one_complete(filesystem, backend):
    """Target runs if ANY parent module is incomplete."""
    shared_target = Target(
        name="SharedTarget",
        inputs=[],
        outputs=["shared.txt"],
        options={},
        working_dir="/some/dir",
    )
    target_a = Target(
        name="TargetA",
        inputs=["shared.txt"],
        outputs=["out_a.txt"],
        options={},
        working_dir="/some/dir",
    )
    target_b = Target(
        name="TargetB",
        inputs=["shared.txt"],
        outputs=["out_b.txt"],
        options={},
        working_dir="/some/dir",
    )

    module_a = Module(name="ModuleA", targets=[shared_target, target_a])
    module_b = Module(name="ModuleB", targets=[shared_target, target_b])

    filesystem.add_file("/some/dir/shared.txt", changed_at=1)
    filesystem.add_file("/some/dir/out_a.txt", changed_at=2)

    graph = Graph.from_targets([module_a, module_b], filesystem)
    target_states = get_status_map(
        graph=graph,
        endpoints=[target_a, target_b],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )

    assert target_states[shared_target] == Status.COMPLETED

    assert target_states[target_a] == Status.COMPLETED

    assert target_states[target_b] == Status.SHOULDRUN


def test_scheduling_target_in_multiple_modules_all_complete(filesystem, backend):
    """Target is skipped if ALL parent modules are complete."""
    shared_target = Target(
        name="SharedTarget",
        inputs=[],
        outputs=["shared.txt"],
        options={},
        working_dir="/some/dir",
    )
    target_a = Target(
        name="TargetA",
        inputs=["shared.txt"],
        outputs=["out_a.txt"],
        options={},
        working_dir="/some/dir",
    )
    target_b = Target(
        name="TargetB",
        inputs=["shared.txt"],
        outputs=["out_b.txt"],
        options={},
        working_dir="/some/dir",
    )

    module_a = Module(name="ModuleA", targets=[shared_target, target_a])
    module_b = Module(name="ModuleB", targets=[shared_target, target_b])

    filesystem.add_file("/some/dir/shared.txt", changed_at=1)
    filesystem.add_file("/some/dir/out_a.txt", changed_at=2)
    filesystem.add_file("/some/dir/out_b.txt", changed_at=2)

    graph = Graph.from_targets([module_a, module_b], filesystem)
    target_states = get_status_map(
        graph=graph,
        endpoints=[target_a, target_b],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
    )

    assert target_states[shared_target] == Status.COMPLETED
    assert target_states[target_a] == Status.COMPLETED
    assert target_states[target_b] == Status.COMPLETED


def test_scheduling_respects_modules_false(filesystem, backend):
    """With respect_modules=False, module completion is ignored."""
    target1 = Target(
        name="Target1",
        inputs=["input.txt"],
        outputs=["out1.txt"],
        options={},
        working_dir="/some/dir",
    )
    module1 = Module(name="Module1", targets=[target1])

    filesystem.add_file("/some/dir/input.txt", changed_at=2)
    filesystem.add_file("/some/dir/out1.txt", changed_at=1)

    graph = Graph.from_targets([module1], filesystem)

    target_states = get_status_map(
        graph=graph,
        endpoints=[target1],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
        respect_modules=True,
    )
    assert target_states[target1] == Status.COMPLETED

    target_states = get_status_map(
        graph=graph,
        endpoints=[target1],
        fs=filesystem,
        backend=backend,
        spec_hashes=NoopSpecHashes(),
        respect_modules=False,
    )
    assert target_states[target1] == Status.SHOULDRUN


def test_override_spec_hashes(backend, spec_hashes, filesystem):
    filesystem.add_file("/some/dir/input.txt", changed_at=1)
    filesystem.add_file("/some/dir/output.txt", changed_at=2)

    target = Target(
        "TestTarget",
        inputs=["input.txt"],
        outputs=["output.txt"],
        options={},
        working_dir="/some/dir",
        spec="foo",
    )
    graph = Graph.from_targets([target], filesystem)

    target_states = get_status_map(
        graph=graph,
        fs=filesystem,
        backend=backend,
        spec_hashes=spec_hashes,
    )
    assert target_states[target] == Status.SHOULDRUN

    spec_hashes.update(target)

    target_states = get_status_map(
        graph=graph,
        fs=filesystem,
        backend=backend,
        spec_hashes=spec_hashes,
    )
    assert target_states[target] == Status.COMPLETED

    target.override_spec_hash = "manual_override"

    target_states = get_status_map(
        graph=graph,
        fs=filesystem,
        backend=backend,
        spec_hashes=spec_hashes,
    )
    assert target_states[target] == Status.SHOULDRUN
