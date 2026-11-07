from gwf import Workflow
from gwf.core import Graph
from gwf.log_storage import (
    clean_logs,
    get_log_paths,
    get_logs_dir,
    init_log_storage,
    prepare_log_storage_for_target,
)


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


def test_log_paths_keep_dotted_target_names(tmp_path):
    names = [f"sample.{i}" for i in range(300)]
    paths = [path for name in names for path in get_log_paths(tmp_path, name)]

    assert len(set(paths)) == len(paths)
    stdout_path, stderr_path = get_log_paths(tmp_path, "sample.1")
    assert stdout_path.name == "sample.1.stdout"
    assert stderr_path.name == "sample.1.stderr"


def test_init_log_storage_migrates_dotted_target_names(tmp_path):
    logs_dir = get_logs_dir(tmp_path)
    logs_dir.mkdir(parents=True)
    logs_dir.joinpath("sample.1.stdout").write_text("out")
    logs_dir.joinpath("sample.1.stderr").write_text("err")

    init_log_storage(tmp_path)

    stdout_path, stderr_path = get_log_paths(tmp_path, "sample.1")
    assert stdout_path.read_text() == "out"
    assert stderr_path.read_text() == "err"
