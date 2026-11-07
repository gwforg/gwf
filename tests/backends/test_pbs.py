import json

import pytest

from gwf import Target
from gwf.backends import pbs
from gwf.backends.base import BackendStatus
from gwf.backends.exceptions import BackendError


@pytest.fixture
def fake_pbs(monkeypatch):
    calls = []

    def fake_call(executable, *args, input=None, environ=None):
        calls.append((executable, args))
        if executable == "qsub":
            return "1234.server\n"
        if executable == "qstat":
            job_id = args[-1]
            if job_id == "1234":
                return "Job Id: 1234.server\n    job_state = R\n"
            if job_id == "999":
                raise BackendError("qstat: Unknown Job Id 999.server\n")
            if job_id == "998":
                raise BackendError(
                    "qstat: 998.server Job has finished, use -x or -H to obtain "
                    "historical job information\n"
                )
            raise BackendError("qstat: cannot connect to server\n")
        raise AssertionError(f"unexpected call to {executable}")

    monkeypatch.setattr(pbs, "call", fake_call)
    return calls


def _target(name):
    return Target(
        name,
        inputs=[],
        outputs=[],
        options=dict(pbs.TARGET_DEFAULTS),
        spec="echo hello",
    )


def _track(working_dir, jobs):
    working_dir.joinpath(".gwf", "logs").mkdir(parents=True)
    working_dir.joinpath(".gwf", "pbs-backend-tracked.json").write_text(
        json.dumps(jobs)
    )


def test_finished_jobs_unknown_to_qstat_have_unknown_status(tmp_path, fake_pbs):
    _track(tmp_path, {"Running": "1234", "Purged": "999", "Finished": "998"})

    with pbs.create_backend(str(tmp_path)) as backend:
        assert backend.status(_target("Running")) == BackendStatus.RUNNING
        assert backend.status(_target("Purged")) == BackendStatus.UNKNOWN
        assert backend.status(_target("Finished")) == BackendStatus.UNKNOWN


def test_other_qstat_errors_are_raised(tmp_path, fake_pbs):
    _track(tmp_path, {"Broken": "777"})

    with pytest.raises(BackendError, match="cannot connect"):
        pbs.create_backend(str(tmp_path))
