from gwf.backends import lsf
from gwf.backends.base import BackendStatus


def test_bjobs_states_are_parsed_and_unknown_states_tolerated(monkeypatch):
    def fake_call(executable, *args, input=None, environ=None):
        assert executable == "bjobs"
        return "1,RUN\n2,PROV\n3,SOME_NEW_STATE\n4,DONE\n"

    monkeypatch.setattr(lsf, "call", fake_call)
    ops = lsf.LSFOps(".", target_defaults=lsf.TARGET_DEFAULTS)

    assert ops.get_job_states(["1", "2", "3"]) == {
        "1": BackendStatus.RUNNING,
        "2": BackendStatus.SUBMITTED,
        "3": BackendStatus.UNKNOWN,
    }
