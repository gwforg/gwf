from unittest import TestCase
from unittest.mock import patch

from click.testing import CliRunner

from exceptions import NoLogFoundError
from gwf.backends import Backend, Status


def touch_file(path, contents=None):
    with open(path, 'w') as fileobj:
        if contents is not None:
            fileobj.write(contents)


class GWFTestCase(TestCase):

    def create_patch(self, name, **kwargs):
        patcher = patch(name, **kwargs)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing


class CliTestCase(TestCase):

    def setUp(self):
        self.runner = CliRunner()

        fs = self.runner.isolated_filesystem()
        fs.__enter__()
        self.addCleanup(fs.__exit__, None, None, None)


class MockBackend(Backend):

    def __init__(self):
        self._tracked = {}

    def submit(self, target, dependencies):
        self._tracked[target] = Status.SUBMITTED

    def cancel(self, target):
        del self._tracked[target]

    def status(self, target):
        return self._tracked.get(target, Status.UNKNOWN)

    def logs(self, target, stderr=False):
        raise NoLogFoundError

    def close(self):
        pass

    def set_status(self, target, status):
        assert status in (Status.RUNNING, Status.UNKNOWN)
        self._tracked[target] = status
