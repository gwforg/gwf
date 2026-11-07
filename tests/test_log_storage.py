from gwf import Workflow
from gwf.core import Graph
from gwf.log_storage import clean_logs, get_log_paths, prepare_log_storage_for_target


def _write_logs(working_dir, name):
    prepare_log_storage_for_target(working_dir, name)
    stdout_path, stderr_path = get_log_paths(working_dir, name)
    stdout_path.write_text("out")
    stderr_path.write_text("err")
    return stdout_path, stderr_path


def test_clean_logs_keeps_active_and_removes_stale_logs(tmp_path):
    workflow = Workflow(working_dir=str(tmp_path))
    workflow.target("Active", inputs=[], outputs=[])
    graph = Graph.from_targets(workflow.targets)

    active_stdout, active_stderr = _write_logs(tmp_path, "Active")
    stale_stdout, stale_stderr = _write_logs(tmp_path, "Stale")

    clean_logs(str(tmp_path), graph)

    assert active_stdout.exists()
    assert active_stderr.exists()
    assert not stale_stdout.exists()
    assert not stale_stderr.exists()
    # "Active" and "Stale" hash into different directories, so the stale
    # target's directory is emptied and should be removed.
    assert stale_stdout.parent != active_stdout.parent
    assert not stale_stdout.parent.exists()
