import pickle
import sys

import pytest

from gwf.utils import ensure_trailing_newline, is_valid_name, load_workflow


def test_ensure_trailing_newline():
    assert ensure_trailing_newline("") == "\n"
    assert ensure_trailing_newline("foo\nbar\n") == "foo\nbar\n"
    assert ensure_trailing_newline("foo\nbar") == "foo\nbar\n"


@pytest.mark.parametrize("name", ["foo", "_foo", "Foo_1", "sample.1"])
def test_is_valid_name_accepts_valid_names(name):
    assert is_valid_name(name)


@pytest.mark.parametrize("name", ["", "1foo", "foo bar", "foo-bar", "foo\n", "\nfoo"])
def test_is_valid_name_rejects_invalid_names(name):
    assert not is_valid_name(name)


@pytest.fixture
def workflow_file(tmp_path):
    """Write a workflow file with a unique module name and unload it after."""
    loaded = []

    def write(name, source):
        path = tmp_path / f"{name}.py"
        path.write_text(source)
        loaded.append(name)
        return path

    yield write
    for name in loaded:
        sys.modules.pop(name, None)


def test_load_workflow_supports_dataclasses_with_postponed_annotations(
    workflow_file,
):
    path = workflow_file(
        "wf_dataclass",
        """from __future__ import annotations
import dataclasses

@dataclasses.dataclass
class Sample:
    name: str

gwf = Sample("a")
""",
    )

    assert load_workflow(path, "gwf").name == "a"


def test_load_workflow_allows_pickling_objects_defined_in_workflow(workflow_file):
    path = workflow_file(
        "wf_executor",
        """from gwf import Workflow

class MyExecutor:
    def get_command(self, spec_path, workflow_root):
        return [spec_path]

gwf = Workflow(working_dir="/some/dir")
gwf.target("A", inputs=[], outputs=[], executor=MyExecutor())
""",
    )

    target = load_workflow(path, "gwf").targets["A"]

    assert type(pickle.loads(pickle.dumps(target)).executor).__name__ == "MyExecutor"


def test_load_workflow_unregisters_module_that_fails_to_load(workflow_file):
    path = workflow_file("wf_broken", "raise RuntimeError('broken')\n")

    with pytest.raises(RuntimeError, match="broken"):
        load_workflow(path, "gwf")

    assert "wf_broken" not in sys.modules
