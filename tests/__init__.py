import unittest
from unittest.mock import patch

from gwf import Target, Workflow


class GWFTestCase(unittest.TestCase):

    def create_patch(self, name):
        patcher = patch(name)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing
