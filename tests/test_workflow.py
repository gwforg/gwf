from unittest.mock import patch

import pytest

from gwf import AnonymousTarget, Workflow
from gwf.exceptions import TypeError, WorkflowError


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


@patch("gwf.workflow.sys._getframe", autospec=True)
@patch("gwf.workflow.inspect.getfile", return_value="/some/path/file.py", autospec=True)
def test_workflow_computes_working_dir_when_not_initialized_with_working_dir(
    inspect_getfile_mock, sys_getframe_mock
):
    workflow = Workflow()

    assert sys_getframe_mock.call_count == 1
    assert inspect_getfile_mock.call_count == 1
    assert workflow.working_dir == "/some/path"


@patch("gwf.workflow._glob", autospec=True)
def test_glob_with_relative_path_searches_relative_to_working_dir(glob_mock):
    workflow = Workflow(working_dir="/some/path")
    workflow.glob("*.fa")
    glob_mock.assert_called_once_with("/some/path/*.fa")


@patch(
    "gwf.workflow._glob",
    return_value=["/other/path/A.fa", "/other/path/B.fa"],
    autospec=True,
)
def test_glob_with_absolute_path_does_not_search_relative_to_working_dir(glob_mock):
    workflow = Workflow(working_dir="/some/path")
    res = workflow.glob("/other/path/*.fa")
    assert res == ["/other/path/A.fa", "/other/path/B.fa"]
    glob_mock.assert_called_once_with("/other/path/*.fa")


@patch("gwf.workflow._iglob", autospec=True)
def test_iglob_with_relative_path_searches_relative_to_working_dir(iglob_mock):
    workflow = Workflow(working_dir="/some/path")
    workflow.iglob("*.fa")
    iglob_mock.assert_called_once_with("/some/path/*.fa")


@patch(
    "gwf.workflow._iglob",
    return_value=["/other/path/A.fa", "/other/path/B.fa"],
    autospec=True,
)
def test_iglob_with_absolute_path_does_not_search_relative_to_working_dir(iglob_mock):
    workflow = Workflow(working_dir="/some/path")
    res = list(workflow.iglob("/other/path/*.fa"))
    assert res == ["/other/path/A.fa", "/other/path/B.fa"]
    iglob_mock.assert_called_once_with("/other/path/*.fa")


@patch("gwf.workflow.subprocess.check_output", autospec=True)
def test_shell_calls_subprocess_with_same_working_dir_as_workflow_in_a_shell(
    mock_check_output,
):
    workflow = Workflow(working_dir="/some/path")
    workflow.shell("echo hello")
    mock_check_output.assert_called_once_with(
        "echo hello", cwd="/some/path", shell=True
    )


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
def mock_template(
    mocker,
):
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
