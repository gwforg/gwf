import io
from unittest.mock import patch

from gwf import Target
from gwf.backends.testing import TestingBackend


def test_backend_logs():
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    backend = TestingBackend()
    with patch.object(
        backend.log_manager, "open_stdout", return_value=io.StringIO("foo")
    ):
        assert backend.logs(target).read() == "foo"

    with patch.object(
        backend.log_manager, "open_stderr", return_value=io.StringIO("bar")
    ):
        assert backend.logs(target, stderr=True).read() == "bar"
