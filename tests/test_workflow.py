from pathlib import Path

import pytest

from gwf import AnonymousTarget, Workflow
from gwf.exceptions import WorkflowError


def test_target_with_no_input_has_empty_inputs_attribute():
    workflow = Workflow()
    target = workflow.target("TestTarget", inputs=[], outputs=[])
    assert target.inputs == []


def test_target_with_no_output_has_empty_outputs_attribute():
    workflow = Workflow()
    target = workflow.target("TestTarget", inputs=[], outputs=[])
    assert target.outputs == []


def test_adding_a_target_makes_it_available_to_the_workflow():
    workflow = Workflow()
    target = workflow.target("TestTarget", inputs=[], outputs=[])

    assert "TestTarget" in workflow.targets
    assert target in workflow.targets.values()


def test_adding_two_targets_with_the_same_names_should_raise_an_exception():
    workflow = Workflow()
    workflow.target("TestTarget", inputs=[], outputs=[])
    with pytest.raises(WorkflowError):
        workflow.target("TestTarget", inputs=[], outputs=[])


def test_target_from_template_returning_anonymous_target():
    def template_returning_anonymous_target():
        return AnonymousTarget(
            inputs=[],
            outputs=[],
            options={},
            working_dir="/some/dir",
            spec="this is the spec",
        )

    workflow = Workflow(working_dir="/some/dir")
    workflow.target_from_template("TestTarget", template_returning_anonymous_target())
    assert "TestTarget" in workflow.targets


def test_target_from_template_returning_anonymous_target_without_working_dir():
    def template_returning_anonymous_target_without_working_dir():
        return AnonymousTarget(
            inputs=["hello.txt"], outputs=[], options={}, spec="this is the spec"
        )

    workflow = Workflow(working_dir="/some/dir")
    workflow.target_from_template(
        "TestTarget", template_returning_anonymous_target_without_working_dir()
    )
    assert "TestTarget" in workflow.targets


def test_targets_inherit_workflow_working_dir_with_given_working_dir():
    workflow = Workflow(working_dir="/some/path")
    target = workflow.target("TestTarget", inputs=[], outputs=[])
    assert target.working_dir == "/some/path"


def test_targets_inherit_workflow_defaults():
    workflow = Workflow(defaults={"cores": 8, "memory": "8g"})
    target = workflow.target("TestTarget", inputs=[], outputs=[])
    assert target.options == {"cores": 8, "memory": "8g"}


def test_target_options_override_defaults():
    workflow = Workflow(defaults={"cores": 8, "memory": "8g"})
    target = workflow.target("TestTarget", inputs=[], outputs=[], cores=16)
    assert target.options == {"cores": 16, "memory": "8g"}


@pytest.mark.parametrize(
    "option_name", ["slurm_args", "sge_args", "lsf_args", "pbs_args"]
)
def test_targets_inherit_workflow_submission_argument_defaults(option_name):
    args = ["--unsupported-option", "value"]
    workflow = Workflow(defaults={option_name: args})

    target = workflow.target("TestTarget", inputs=[], outputs=[])

    assert target.options == {option_name: args}


@pytest.mark.parametrize(
    "option_name", ["slurm_args", "sge_args", "lsf_args", "pbs_args"]
)
def test_target_submission_arguments_override_workflow_defaults(option_name):
    workflow = Workflow(defaults={option_name: ["--default-option"]})

    target = workflow.target(
        "TestTarget",
        inputs=[],
        outputs=[],
        **{option_name: ["--target-option", "value"]},
    )

    assert target.options == {option_name: ["--target-option", "value"]}


def test_workflow_computes_working_dir_when_not_initialized_with_working_dir():
    workflow = Workflow()

    assert workflow.working_dir == str(Path(__file__).parent)


@pytest.mark.parametrize("method", ["glob", "iglob"])
def test_glob_with_relative_path_searches_relative_to_working_dir(tmp_path, method):
    expected = tmp_path / "sequence.fa"
    expected.touch()
    (tmp_path / "sequence.txt").touch()

    workflow = Workflow(working_dir=tmp_path)

    assert list(getattr(workflow, method)("*.fa")) == [str(expected)]


@pytest.mark.parametrize("method", ["glob", "iglob"])
def test_glob_with_absolute_path_does_not_search_relative_to_working_dir(
    tmp_path, method
):
    working_dir = tmp_path / "workflow"
    other_dir = tmp_path / "inputs"
    working_dir.mkdir()
    other_dir.mkdir()
    expected = other_dir / "sequence.fa"
    expected.touch()

    workflow = Workflow(working_dir=working_dir)

    assert list(getattr(workflow, method)(str(other_dir / "*.fa"))) == [str(expected)]


def test_shell_runs_command_in_workflow_working_dir(tmp_path):
    workflow = Workflow(working_dir=tmp_path)

    assert workflow.shell("printf hello > greeting.txt") == b""
    assert (tmp_path / "greeting.txt").read_text() == "hello"


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


def test_map_passes_string_inputs_to_template():
    files = ["a", "b", "c"]

    def template(path):
        return AnonymousTarget(inputs={"path": path}, outputs=[], options={})

    workflow = Workflow(working_dir="/some/dir")
    targets = workflow.map(template, files)

    assert [target.inputs for target in targets] == [
        {"path": "a"},
        {"path": "b"},
        {"path": "c"},
    ]


def test_map_unpacks_tuple_inputs_for_template():
    files = [("a", "/foo"), ("b", "/foo"), ("c", "/foo")]

    def template(path, output_dir):
        return AnonymousTarget(
            inputs={"path": path, "output_dir": output_dir}, outputs=[], options={}
        )

    workflow = Workflow(working_dir="/some/dir")
    targets = workflow.map(template, files)

    assert [target.inputs for target in targets] == [
        {"path": "a", "output_dir": "/foo"},
        {"path": "b", "output_dir": "/foo"},
        {"path": "c", "output_dir": "/foo"},
    ]


def test_map_unpacks_mapping_inputs_for_template():
    files = [
        {"path": "a", "output_dir": "foo/"},
        {"path": "b", "output_dir": "foo/"},
        {"path": "c", "output_dir": "foo/"},
    ]

    def template(path, output_dir):
        return AnonymousTarget(
            inputs={"path": path, "output_dir": output_dir}, outputs=[], options={}
        )

    workflow = Workflow(working_dir="/some/dir")
    targets = workflow.map(template, files)

    assert [target.inputs for target in targets] == files


def test_map_passes_extra_arguments_to_template():
    files = [{"path": "a"}, {"path": "b"}, {"path": "c"}]

    def template(path, output_dir):
        return AnonymousTarget(
            inputs={"path": path, "output_dir": output_dir}, outputs=[], options={}
        )

    workflow = Workflow(working_dir="/some/dir")
    targets = workflow.map(template, files, extra={"output_dir": "foo/"})

    assert [target.inputs for target in targets] == [
        {"path": "a", "output_dir": "foo/"},
        {"path": "b", "output_dir": "foo/"},
        {"path": "c", "output_dir": "foo/"},
    ]


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
