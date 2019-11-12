import pytest
from unittest.mock import Mock, patch

from gwf.workflow import AnonymousTarget, Workflow
from gwf.exceptions import NameError, TypeError, WorkflowError


def test_workflow_with_invalid_name_raises_error():
    with pytest.raises(NameError):
        Workflow(name="123abc")


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
