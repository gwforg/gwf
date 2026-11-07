from gwf.backends import slurm
from gwf.backends.base import BackendStatus


def test_sacct_states_are_parsed_and_unknown_states_tolerated(monkeypatch):
    def fake_call(executable, *args, input=None, environ=None):
        assert executable == "sacct"
        return "1|RUNNING\n2|CANCELLED by 1000\n3|SOME_NEW_STATE\n"

    monkeypatch.setattr(slurm, "call", fake_call)
    ops = slurm.SlurmOps(".", "full", True, target_defaults=slurm.TARGET_DEFAULTS)

    assert ops.get_job_states_from_sacct(["1", "2", "3"]) == {
        "1": BackendStatus.RUNNING,
        "2": BackendStatus.CANCELLED,
        "3": BackendStatus.UNKNOWN,
    }
