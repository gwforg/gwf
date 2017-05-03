from unittest import TestCase

from gwf import Target
from gwf.backends.testing import TestingBackend


class TestBackendType(TestCase):

    def setUp(self):
        self.backend = TestingBackend(working_dir='/some/dir')

    def test_inject_target_defaults_into_target_options_on_submit(self):
        target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        self.backend.submit(target, dependencies=[])
        self.assertEqual(target.options, {'cores': 2, 'memory': '18g'})

        target = Target('TestTarget', inputs=[], outputs=[], options={'cores': 32}, working_dir='/some/dir')
        self.backend.submit(target, dependencies=[])
        self.assertEqual(target.options, {'cores': 32, 'memory': '18g'})

    def test_warn_user_when_submitting_target_with_unsupported_option(self):
        target = Target('TestTarget', inputs=[], outputs=[], options={'foo': 'bar'}, working_dir='/some/dir')
        with self.assertWarns(Warning):
            self.backend.submit(target, dependencies=[])

