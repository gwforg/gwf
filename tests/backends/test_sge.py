import pytest

from gwf import Target
from gwf.backends import sge
from gwf.backends.base import BackendStatus

QSTAT_XML = """<?xml version='1.0'?>
<job_info>
  <queue_info>
    <job_list state="running">
      <JB_job_number>{job_id}</JB_job_number>
      <state>r</state>
    </job_list>
  </queue_info>
</job_info>
"""


@pytest.fixture
def fake_sge(monkeypatch):
    calls = []
    job_ids = iter(["1234", "1235"])

    def fake_call(executable, *args, input=None, environ=None):
        calls.append((executable, args))
        if executable == "qsub":
            # Real `qsub -terse` output ends with a newline.
            return f"{next(job_ids)}\n"
        if executable == "qstat":
            return QSTAT_XML.format(job_id="1234")
        raise AssertionError(f"unexpected call to {executable}")

    monkeypatch.setattr(sge, "call", fake_call)
    return calls


def _target(name):
    return Target(
        name,
        inputs=[],
        outputs=[],
        options=dict(sge.TARGET_DEFAULTS),
        spec="echo hello",
    )


def test_submitted_job_id_matches_qstat_job_id(tmp_path, fake_sge):
    tmp_path.joinpath(".gwf").mkdir()
    target = _target("Foo")

    with sge.create_backend(str(tmp_path)) as backend:
        backend.submit(target, dependencies=[])
        assert backend.get_tracked_id(target) == "1234"

    with sge.create_backend(str(tmp_path)) as backend:
        assert backend.status(target) == BackendStatus.RUNNING


def test_dependency_ids_are_passed_without_newlines(tmp_path, fake_sge):
    tmp_path.joinpath(".gwf").mkdir()
    upstream, downstream = _target("Upstream"), _target("Downstream")

    with sge.create_backend(str(tmp_path)) as backend:
        backend.submit(upstream, dependencies=[])
        backend.submit(downstream, dependencies=[upstream])

    _, args = [call for call in fake_sge if call[0] == "qsub"][-1]
    assert args[args.index("-hold_jid") + 1] == "1234"
