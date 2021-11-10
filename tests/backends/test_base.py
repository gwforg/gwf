import io
import logging
from unittest.mock import patch

from gwf import Target
from gwf.backends.testing import TestingBackend


def test_backend_submit_full_injects_backend_defaults(backend):
    target1 = Target(
        "TestTarget1", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    target2 = Target(
        "TestTarget2",
        inputs=[],
        outputs=[],
        options={"cores": 32},
        working_dir="/some/dir",
    )

    backend.submit_full(target1, set())
    assert target1.options == {"cores": 1, "memory": "1g"}

    backend.submit_full(target2, set())
    assert target2.options == {"cores": 32, "memory": "1g"}


def test_backend_submit_full_warns_user_when_submitting_target_with_unsupported_option(
    backend, caplog
):
    target = Target(
        "TestTarget",
        inputs=[],
        outputs=[],
        options={"foo": "bar"},
        working_dir="/some/dir",
    )

    backend.submit_full(target, set())

    assert target.options == {"cores": 1, "memory": "1g"}
    assert caplog.record_tuples == [
        (
            "gwf.backends.base",
            logging.WARNING,
            "Option 'foo' used in 'TestTarget' is not supported by backend. Ignored.",
        )
    ]


def test_backend_submit_full_removes_options_with_none_value(backend):
    target = Target(
        "TestTarget",
        inputs=[],
        outputs=[],
        options={"cores": None},
        working_dir="/some/dir",
    )

    backend.submit_full(target, set())
    assert target.options == {"memory": "1g"}


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
