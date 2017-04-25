from unittest import TestCase
from unittest.mock import patch


class GWFTestCase(TestCase):

    def create_patch(self, name):
        patcher = patch(name)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing
