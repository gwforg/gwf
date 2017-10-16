import io
import logging
import os.path
from unittest.mock import patch

import pytest

from gwf import Target
from gwf.backends import Backend
from gwf.backends.testing import TestingBackend
from gwf.exceptions import BackendError
from gwf.utils import PersistableDict


@pytest.fixture
def backend():
    return TestingBackend()


def test_inject_target_defaults_into_target_options_on_submit(backend):
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')
    backend.submit(target, dependencies=[])
    assert target.options == {'cores': 2, 'memory': '18g'}

    target = Target('TestTarget', inputs=[], outputs=[], options={'cores': 32}, working_dir='/some/dir')
    backend.submit(target, dependencies=[])
    assert target.options == {'cores': 32, 'memory': '18g'}


def test_warn_user_when_submitting_target_with_unsupported_option(backend, caplog):
    target = Target('TestTarget', inputs=[], outputs=[], options={'foo': 'bar'}, working_dir='/some/dir')
    backend.submit(target, dependencies=[])
    assert caplog.record_tuples == [
        ('gwf.backends.base', logging.WARNING, 'Option "foo" used in "TestTarget" is not supported by backend. Ignored.'),
    ]


def test_remove_options_with_none_value(backend):
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')
    backend.submit(target, dependencies=[])
    assert target.options == {'cores': 2, 'memory': '18g'}


def test_backend_logs():
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

    backend = TestingBackend()
    with patch.object(backend.log_manager, 'open_stdout', return_value=io.StringIO('foo')):
        assert backend.logs(target).read() == 'foo'

    with patch.object(backend.log_manager, 'open_stderr', return_value=io.StringIO('bar')):
        assert backend.logs(target, stderr=True).read() == 'bar'


def test_raise_exception_if_backend_does_not_implement_all_methods():
    with pytest.raises(BackendError):
        class InvalidTestingBackend(Backend):
            def submit(self, target, dependencies):
                pass

            def cancel(self, target):
                pass

            def status(self, target):
                pass


def test_dict_is_empty_if_file_does_not_exist():
    d = PersistableDict('test.json')
    assert d == {}


def test_dict_write_read():
    d1 = PersistableDict('/tmp/test.json')
    d1['foo'] = 'bar'
    d1.persist()

    assert os.path.exists('/tmp/test.json')

    d2 = PersistableDict('/tmp/test.json')
    assert d2 == {'foo': 'bar'}

