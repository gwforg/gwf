"""Commands that change nothing, like `gwf status`, must not write state files."""

import json

from gwf import Target
from gwf.backends.base import TrackingBackend
from gwf.core import FileSpecHashes


class FakeOps:
    def __init__(self):
        self.target_defaults = {}

    def get_job_states(self, tracked_jobs):
        return {}

    def submit_target(self, target, dependencies):
        return f"{target.name}-1"

    def close(self):
        pass


def _target(name, spec=""):
    return Target(name, inputs=[], outputs=[], options={}, spec=spec)


def _read(path):
    return json.loads(path.read_text())


def test_tracking_backend_without_submissions_does_not_overwrite_state(tmp_path):
    (tmp_path / ".gwf").mkdir()

    # `gwf status` and `gwf run` start at the same time.
    status_backend = TrackingBackend(str(tmp_path), name="fake", ops=FakeOps())
    with TrackingBackend(str(tmp_path), name="fake", ops=FakeOps()) as run_backend:
        run_backend.submit(_target("A"), dependencies=[])

    # `gwf status` finishes last.
    status_backend.close()

    assert _read(tmp_path / ".gwf" / "fake-backend-tracked.json") == {"A": "A-1"}


def test_spec_hashes_without_changes_do_not_overwrite_file(tmp_path):
    path = tmp_path / "spec-hashes.json"

    status_hashes = FileSpecHashes(str(path))
    with FileSpecHashes(str(path)) as run_hashes:
        run_hashes.update(_target("A", spec="echo a"))
    status_hashes.close()

    assert set(_read(path)) == {"A"}


def test_spec_hashes_invalidate_is_persisted(tmp_path):
    path = tmp_path / "spec-hashes.json"
    with FileSpecHashes(str(path)) as hashes:
        hashes.update(_target("A", spec="echo a"))
        hashes.update(_target("B", spec="echo b"))

    with FileSpecHashes(str(path)) as hashes:
        hashes.invalidate(_target("A"))

    assert set(_read(path)) == {"B"}
