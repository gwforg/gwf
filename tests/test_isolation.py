import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

import pytest

from gwf import AnonymousTarget, IsolationConfig, Workflow
from gwf.core import Target
from gwf.exec.isolation import IsolationError
from gwf.executors import render_exec_script


def run_target(tmp_path, *, inputs, outputs, isolation, spec):
    target = Target(
        "IsolatedTarget",
        inputs=inputs,
        outputs=outputs,
        options={},
        isolation=isolation,
        working_dir=str(tmp_path),
        spec=spec,
    )
    with render_exec_script(target) as script:
        pass

    script_path = tmp_path / "job-script"
    script_path.write_text(script.getvalue())
    script_path.chmod(0o700)

    environ = os.environ.copy()
    environ["GWF_EXEC_WORKFLOW_ROOT"] = str(tmp_path)
    return subprocess.run(
        [script_path],
        cwd=tmp_path,
        env=environ,
        text=True,
        capture_output=True,
        check=False,
    )


def network_namespaces_available():
    if sys.platform != "linux":
        return False
    try:
        result = subprocess.run(
            ["unshare", "--user", "--map-current-user", "--net", "--", "true"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except FileNotFoundError:
        return False
    return result.returncode == 0


def test_target_inherits_isolation_from_workflow_defaults(tmp_path):
    scratch = tmp_path / "scratch"
    workflow = Workflow(
        working_dir=str(tmp_path),
        defaults={
            "cores": 2,
            "isolation": IsolationConfig(
                inputs="symlink", outputs="copy", root=scratch, network=True
            ),
        },
    )

    target = workflow.target("Target", inputs=[], outputs=[])

    assert isinstance(target.isolation, IsolationConfig)
    assert target.isolation.inputs == "symlink"
    assert target.isolation.outputs == "copy"
    assert target.isolation.root == scratch
    assert target.isolation.network
    assert target.options == {"cores": 2}


def test_target_isolation_overrides_workflow_default(tmp_path):
    workflow = Workflow(
        working_dir=str(tmp_path),
        defaults={"isolation": IsolationConfig(inputs="symlink", outputs="copy")},
    )
    isolation = IsolationConfig()

    target = workflow.target("Target", inputs=[], outputs=[], isolation=isolation)

    assert target.isolation is isolation
    assert target.isolation.inputs == "symlink"
    assert target.isolation.outputs == "move"
    assert "isolation" not in target.options


def test_target_isolation_defaults_to_none(tmp_path):
    workflow = Workflow(working_dir=str(tmp_path))

    target = workflow.target("Target", inputs=[], outputs=[])

    assert target.isolation is None


def test_true_enables_default_isolation(tmp_path):
    workflow = Workflow(working_dir=str(tmp_path))

    target = workflow.target("Target", inputs=[], outputs=[], isolation=True)

    assert target.isolation == IsolationConfig()


def test_true_enables_default_isolation_in_workflow_defaults(tmp_path):
    workflow = Workflow(
        working_dir=str(tmp_path),
        defaults={"isolation": True},
    )

    target = workflow.target("Target", inputs=[], outputs=[])

    assert target.isolation == IsolationConfig()


def test_template_isolation_overrides_workflow_default(tmp_path):
    workflow = Workflow(
        working_dir=str(tmp_path),
        defaults={"isolation": IsolationConfig(inputs="symlink", outputs="copy")},
    )
    template = AnonymousTarget(
        inputs=[],
        outputs=[],
        options={"cores": 2},
        isolation=IsolationConfig(inputs="copy", outputs="move"),
    )

    target = workflow.target_from_template("Target", template)

    assert isinstance(target.isolation, IsolationConfig)
    assert target.isolation.inputs == "copy"
    assert target.isolation.outputs == "move"
    assert target.options == {"cores": 2}


def test_template_can_configure_isolation_in_its_options(tmp_path):
    workflow = Workflow(working_dir=str(tmp_path))
    template = AnonymousTarget(
        inputs=[],
        outputs=[],
        options={"cores": 2, "isolation": IsolationConfig(inputs="symlink")},
    )

    target = workflow.target_from_template("Target", template)

    assert isinstance(target.isolation, IsolationConfig)
    assert target.isolation.inputs == "symlink"
    assert target.options == {"cores": 2}


def test_target_isolation_overrides_template_isolation(tmp_path):
    workflow = Workflow(
        working_dir=str(tmp_path), defaults={"isolation": IsolationConfig()}
    )
    template = AnonymousTarget(
        inputs=[], outputs=[], options={}, isolation=IsolationConfig(inputs="symlink")
    )
    isolation = IsolationConfig(outputs="copy")

    target = workflow.target_from_template("Target", template, isolation=isolation)

    assert target.isolation is isolation
    assert "isolation" not in target.options


@pytest.mark.parametrize(
    "isolation",
    [
        False,
        {"inputs": "copy"},
        "copy",
    ],
)
def test_isolation_requires_an_isolation_config(tmp_path, isolation):
    workflow = Workflow(working_dir=str(tmp_path))

    with pytest.raises(TypeError):
        workflow.target("Target", inputs=[], outputs=[], isolation=isolation)


@pytest.mark.parametrize(
    "config",
    [
        lambda: IsolationConfig(inputs="hardlink"),
        lambda: IsolationConfig(outputs="delete"),
    ],
)
def test_invalid_isolation_config_mode_is_rejected(config):
    with pytest.raises(IsolationError):
        config()


def test_isolation_root_must_be_a_path():
    with pytest.raises(TypeError):
        IsolationConfig(root=42)


def test_true_enables_default_isolation_on_anonymous_target():
    target = AnonymousTarget(inputs=[], outputs=[], options={}, isolation=True)

    assert target.isolation == IsolationConfig()


def test_isolation_blocks_network_by_default():
    assert not IsolationConfig().network


def test_isolation_can_allow_network_access():
    assert IsolationConfig(network=True).network


@pytest.mark.skipif(
    not network_namespaces_available(),
    reason="unprivileged network namespaces are unavailable",
)
def test_isolation_runs_target_in_a_private_network_namespace(tmp_path):
    host_namespace = os.readlink("/proc/self/ns/net")

    result = run_target(
        tmp_path,
        inputs=[],
        outputs=["network-namespace.txt"],
        isolation=IsolationConfig(),
        spec="readlink /proc/self/ns/net > network-namespace.txt",
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "network-namespace.txt").read_text().strip() != host_namespace


def test_disabled_isolation_uses_original_working_directory(tmp_path):
    (tmp_path / "undeclared.txt").write_text("visible\n")

    result = run_target(
        tmp_path,
        inputs=[],
        outputs=["output.txt"],
        isolation=None,
        spec="test -e undeclared.txt; pwd > output.txt",
    )

    assert result.returncode == 0, result.stderr
    assert Path((tmp_path / "output.txt").read_text().strip()).resolve() == (
        tmp_path.resolve()
    )


@pytest.mark.parametrize("input_mode", ["copy", "symlink"])
@pytest.mark.parametrize("output_mode", ["copy", "move"])
def test_isolated_target_stages_files_and_runs_in_clean_directory(
    tmp_path, input_mode, output_mode
):
    input_path = tmp_path / "inputs" / "source.txt"
    input_path.parent.mkdir()
    input_path.write_text("input data\n")
    (tmp_path / "undeclared.txt").write_text("must not be visible\n")
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "output.txt").write_text("previous result\n")

    link_assertion = "-L" if input_mode == "symlink" else "! -L"
    result = run_target(
        tmp_path,
        inputs={"source": "inputs/source.txt"},
        outputs={
            "result": "results/output.txt",
            "working_dir": "results/cwd.txt",
        },
        isolation=IsolationConfig(inputs=input_mode, outputs=output_mode, network=True),
        spec=f"""
test ! -e undeclared.txt
test {link_assertion} inputs/source.txt
cat inputs/source.txt > results/output.txt
pwd > results/cwd.txt
""",
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "results" / "output.txt").read_text() == "input data\n"

    isolated_directory = Path((tmp_path / "results" / "cwd.txt").read_text().strip())
    assert isolated_directory.resolve() != tmp_path.resolve()
    assert not isolated_directory.exists()


@pytest.mark.parametrize("configured_root", ["scratch", "absolute"])
def test_isolation_uses_configured_temporary_root(tmp_path, configured_root):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    root = scratch if configured_root == "absolute" else configured_root

    result = run_target(
        tmp_path,
        inputs=[],
        outputs=["working-directory.txt"],
        isolation=IsolationConfig(root=root, network=True),
        spec="pwd > working-directory.txt",
    )

    assert result.returncode == 0, result.stderr
    isolated_working_directory = Path(
        (tmp_path / "working-directory.txt").read_text().strip()
    )
    assert isolated_working_directory.parent.parent == scratch.resolve()
    assert not isolated_working_directory.exists()


def test_missing_isolation_root_is_reported(tmp_path):
    missing_root = tmp_path / "missing"

    result = run_target(
        tmp_path,
        inputs=[],
        outputs=[],
        isolation=IsolationConfig(root=missing_root, network=True),
        spec="true",
    )

    assert result.returncode != 0
    assert "Could not create an isolated working directory" in result.stderr
    assert str(missing_root) in result.stderr


def test_copied_input_cannot_modify_original_file(tmp_path):
    input_path = tmp_path / "input.txt"
    input_path.write_text("original\n")

    result = run_target(
        tmp_path,
        inputs=["input.txt"],
        outputs=["output.txt"],
        isolation=IsolationConfig(inputs="copy", network=True),
        spec="printf 'changed\\n' > input.txt; cp input.txt output.txt",
    )

    assert result.returncode == 0, result.stderr
    assert input_path.read_text() == "original\n"
    assert (tmp_path / "output.txt").read_text() == "changed\n"


def test_failed_target_does_not_publish_partial_outputs(tmp_path):
    output_path = tmp_path / "output.txt"
    output_path.write_text("previous result\n")

    result = run_target(
        tmp_path,
        inputs=[],
        outputs=["output.txt"],
        isolation=IsolationConfig(network=True),
        spec="printf 'partial result\\n' > output.txt; exit 13",
    )

    assert result.returncode == 13
    assert output_path.read_text() == "previous result\n"


def test_transfer_failure_can_leave_partial_outputs(tmp_path):
    output_path = tmp_path / "first.txt"
    output_path.write_text("previous result\n")
    (tmp_path / "second.txt").mkdir()

    result = run_target(
        tmp_path,
        inputs=[],
        outputs=["first.txt", "second.txt"],
        isolation=IsolationConfig(network=True),
        spec="""
printf 'new first\\n' > first.txt
printf 'new second\\n' > second.txt
""",
    )

    assert result.returncode != 0
    assert output_path.read_text() == "new first\n"
    assert (tmp_path / "second.txt").is_dir()


def test_successful_target_fails_when_declared_output_is_missing(tmp_path):
    output_path = tmp_path / "complete.txt"
    output_path.write_text("previous result\n")
    result = run_target(
        tmp_path,
        inputs=[],
        outputs=["complete.txt", "missing.txt"],
        isolation=IsolationConfig(network=True),
        spec="printf 'new result\\n' > complete.txt",
    )

    assert result.returncode != 0
    assert "missing.txt" in result.stderr
    assert output_path.read_text() == "previous result\n"


def test_output_must_remain_inside_isolated_working_directory(tmp_path):
    result = run_target(
        tmp_path,
        inputs=[],
        outputs=["results/output.txt"],
        isolation=IsolationConfig(network=True),
        spec="""
rmdir results
mkdir ../outside
ln -s ../outside results
printf 'escaped\\n' > results/output.txt
""",
    )

    assert result.returncode != 0
    assert "outside its isolated working directory" in result.stderr
    assert not (tmp_path / "results" / "output.txt").exists()


def test_isolated_input_must_be_a_regular_file(tmp_path):
    (tmp_path / "input").mkdir()

    result = run_target(
        tmp_path,
        inputs=["input"],
        outputs=[],
        isolation=IsolationConfig(network=True),
        spec="true",
    )

    assert result.returncode != 0
    assert "not a file" in result.stderr


def test_isolated_output_must_be_a_regular_file(tmp_path):
    result = run_target(
        tmp_path,
        inputs=[],
        outputs=["output"],
        isolation=IsolationConfig(network=True),
        spec="mkdir output",
    )

    assert result.returncode != 0
    assert "not a regular file" in result.stderr


def test_termination_cleans_isolated_directory_and_stops_target(tmp_path):
    target = Target(
        "IsolatedTarget",
        inputs=[],
        outputs=["output.txt"],
        options={},
        isolation=IsolationConfig(network=True),
        working_dir=str(tmp_path),
        spec="""
sleep 30 &
echo TARGET_STARTED CHILD_PID=$! >&2
wait
touch output.txt
""",
    )
    with render_exec_script(target) as script:
        pass

    script_path = tmp_path / "job-script"
    script_path.write_text(script.getvalue())
    script_path.chmod(0o700)
    environ = os.environ.copy()
    environ["GWF_EXEC_WORKFLOW_ROOT"] = str(tmp_path)
    process = subprocess.Popen(
        [script_path],
        cwd=tmp_path,
        env=environ,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    isolated_directory = None
    child_pid = None
    try:
        for line in process.stderr:
            match = re.search(r"created isolated working directory (.+)", line)
            if match:
                isolated_directory = Path(match.group(1))
            match = re.search(r"TARGET_STARTED CHILD_PID=(\d+)", line)
            if match:
                child_pid = int(match.group(1))
                break
        else:
            pytest.fail("target exited before it started")

        process.terminate()
        assert process.wait(timeout=10) == 128 + signal.SIGTERM
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()

    assert isolated_directory is not None
    assert not isolated_directory.exists()
    assert not (tmp_path / "output.txt").exists()
    assert child_pid is not None

    deadline = time.monotonic() + 5
    while True:
        try:
            os.kill(child_pid, 0)
        except ProcessLookupError:
            break
        if time.monotonic() >= deadline:
            os.kill(child_pid, signal.SIGKILL)
            pytest.fail(f"target child process {child_pid} was not terminated")
        time.sleep(0.01)


@pytest.mark.parametrize("declared_path", ["../outside.txt", "/outside.txt"])
def test_isolated_target_rejects_paths_outside_working_directory(
    tmp_path, declared_path
):
    result = run_target(
        tmp_path,
        inputs=[],
        outputs=[declared_path],
        isolation=IsolationConfig(network=True),
        spec="true",
    )

    assert result.returncode != 0
    assert declared_path in result.stderr
