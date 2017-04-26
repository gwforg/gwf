from unittest import TestCase
from unittest.mock import patch

from click.testing import CliRunner


def touch_file(path, contents=None):
    with open(path, 'w') as fileobj:
        if contents is not None:
            fileobj.write(contents)


class GWFTestCase(TestCase):

    def create_patch(self, name):
        patcher = patch(name)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing


class CliTestCase(TestCase):

    def setUp(self):
        self.runner = CliRunner()

        fs = self.runner.isolated_filesystem()
        fs.__enter__()
        self.addCleanup(fs.__exit__, None, None, None)
