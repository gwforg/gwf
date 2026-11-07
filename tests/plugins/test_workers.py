import os

from gwf.cli import main
from gwf.plugins import workers


def test_workers_uses_cores_available_to_the_process(
    cli_runner, simple_workflow, monkeypatch
):
    started = {}

    def fake_start_cluster(working_dir, max_cores, max_memory=None):
        started["max_cores"] = max_cores

    monkeypatch.setattr(workers, "start_cluster", fake_start_cluster)
    monkeypatch.setattr(os, "process_cpu_count", lambda: 3, raising=False)

    result = cli_runner.invoke(main, ["-b", "local", "workers"])

    assert result.exit_code == 0, result.output
    assert started["max_cores"] == 3
