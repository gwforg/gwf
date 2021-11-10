from unittest.mock import Mock, patch

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


def test_target_from_template_returning_tuple():
    def template_returning_tuple():
        return [], [], {}, "this is the spec"

    workflow = Workflow(working_dir="/some/dir")
    workflow.target_from_template("TestTarget", template_returning_tuple())
    assert "TestTarget" in workflow.targets


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


def test_target_from_invalid_template():
    def invalid_template():
        return [], []

    workflow = Workflow()

    with pytest.raises(TypeError):
        workflow.target_from_template("TestTarget", 50)

    with pytest.raises(TypeError):
        workflow.target_from_template("TestTarget", invalid_template())


def test_including_workflow_with_no_name_raises_an_exception():
    workflow = Workflow()
    other_workflow = Workflow()
    with pytest.raises(WorkflowError):
        workflow.include(other_workflow)


def test_including_workflow_object_should_extend_including_workflow():
    workflow = Workflow()
    workflow.target("TestTarget1", inputs=[], outputs=[])

    other_workflow = Workflow(name="foo")
    other_workflow.target("TestTarget2", inputs=[], outputs=[])
    other_workflow.target("TestTarget3", inputs=[], outputs=[])

    workflow.include_workflow(other_workflow)

    assert "TestTarget1" in workflow.targets
    assert "foo.TestTarget2" in workflow.targets
    assert "foo.TestTarget3" in workflow.targets


def test_include_target_from_workflow_in_two_different_workflows_():
    w1 = Workflow()
    target = w1.target("MyTarget", inputs=[], outputs=[])

    w3 = Workflow()
    w3.include(w1, namespace="bar")

    w2 = Workflow()
    w2.include(w1, namespace="foo")

    assert target.name == "MyTarget"
    assert "bar.MyTarget" in w3.targets
    assert w3.targets["bar.MyTarget"].name == "bar.MyTarget"
    assert "foo.MyTarget" in w2.targets
    assert w2.targets["foo.MyTarget"].name == "foo.MyTarget"


def test_including_workflow_with_same_name_as_this_workflow_raises_an_exception():
    workflow = Workflow(name="foo")
    other_workflow = Workflow(name="foo")
    with pytest.raises(WorkflowError):
        workflow.include(other_workflow)


@patch("gwf.workflow.load_workflow", autospec=True)
def test_including_workflow_from_path(mock_load_workflow):
    workflow = Workflow()
    workflow.target("TestTarget1", inputs=[], outputs=[])

    other_workflow = Workflow()
    other_workflow.target("TestTarget2", inputs=[], outputs=[])
    other_workflow.target("TestTarget3", inputs=[], outputs=[])

    mock_load_workflow.return_value = other_workflow

    workflow.include_path("/path/to/other_workflow.py", namespace="other")
    assert workflow.targets.keys() == {
        "TestTarget1",
        "other.TestTarget2",
        "other.TestTarget3",
    }


def test_including_workflow_instance_dispatches_to_include_workflow():
    workflow = Workflow()
    other_workflow = Workflow()

    with patch.object(
        workflow, "include_workflow", autospec=True
    ) as mock_include_workflow:
        workflow.include(other_workflow)
        mock_include_workflow.assert_called_once_with(other_workflow, namespace=None)


def test_including_workflow_path_dispatches_to_include_path():
    workflow = Workflow()

    with patch.object(workflow, "include_path", autospec=True) as mock_include_path:
        workflow.include("/path/to/other_workflow.py")
        mock_include_path.assert_called_once_with(
            "/path/to/other_workflow.py", namespace=None
        )


@patch("gwf.workflow.inspect.ismodule", return_value=True, autospec=True)
def test_including_workflow_module_gets_workflow_attribute_and_dispatches_to_include_workflow(
    mock_ismodule,
):
    workflow = Workflow(working_dir="/some/dir")
    other_workflow = Workflow(working_dir="/some/other/dir")

    mock_module = Mock()
    mock_module.gwf = other_workflow

    with patch.object(
        workflow, "include_workflow", autospec=True
    ) as mock_include_workflow:
        workflow.include(mock_module)

        mock_ismodule.assert_called_once_with(mock_module)
        mock_include_workflow.assert_called_once_with(other_workflow, namespace=None)


def test_including_non_module_str_and_object_value_raises_type_error():
    workflow = Workflow(working_dir="/some/dir")
    with pytest.raises(TypeError):
        workflow.include(42)


def test_targets_inherit_workflow_working_dir_with_given_working_dir():
    workflow = Workflow(working_dir="/some/path")
    target = workflow.target("TestTarget", inputs=[], outputs=[])
    target.working_dir == "/some/path"


def test_targets_inherit_workflow_defaults():
    workflow = Workflow(defaults={"cores": 8, "memory": "8g"})
    target = workflow.target("TestTarget", inputs=[], outputs=[])
    target.options == {"cores": 8, "memory": "8g"}


def test_target_options_override_defaults():
    workflow = Workflow(defaults={"cores": 8, "memory": "8g"})
    target = workflow.target("TestTarget", inputs=[], outputs=[], cores=16)
    target.options == {"cores": 16, "memory": "8g"}


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
