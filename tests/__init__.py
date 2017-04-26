from unittest import TestCase
from unittest.mock import patch

from click.testing import CliRunner


class GWFTestCase(TestCase):

    def create_patch(self, name):
        patcher = patch(name)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing


class CliTestCase(TestCase):

    def setUp(self):
        self.runner = CliRunner()
