from unittest.mock import Mock

from gwf import Target
from gwf.backends import Status
from gwf.filtering import StatusFilter, Criteria, filter, NameFilter, EndpointFilter

from tests import GWFTestCase


class TestStatusFilter(GWFTestCase):

    def setUp(self):
        self.target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        self.backend = Mock(name='Backend', spec_set=['status'])
        self.graph = Mock(name='Graph', spec_set=['targets', 'endpoints', 'should_run'])

    def test_should_not_use_filter(self):
        criteria = Criteria(status=None)
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertFalse(f.use())

    def test_should_use_filter(self):
        criteria = Criteria(status='running')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.use())

    def test_status_is_completed_with_completed_target(self):
        self.backend.status.return_value = Status.UNKNOWN
        self.graph.should_run.return_value = False

        criteria = Criteria(status='completed')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.predicate(self.target))

    def test_status_is_completed_with_running_target(self):
        self.backend.status.return_value = Status.RUNNING
        self.graph.should_run.return_value = True

        criteria = Criteria(status='completed')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertFalse(f.predicate(self.target))

    def test_status_is_shouldrun_with_target_that_should_run(self):
        self.backend.status.return_value = Status.UNKNOWN
        self.graph.should_run.return_value = True

        criteria = Criteria(status='shouldrun')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.predicate(self.target))

    def test_status_is_shouldrun_with_target_that_should_not_run(self):
        self.backend.status.return_value = Status.RUNNING
        self.graph.should_run.return_value = False

        criteria = Criteria(status='shouldrun')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertFalse(f.predicate(self.target))

    def test_status_is_running_with_target_that_is_running(self):
        self.backend.status.return_value = Status.RUNNING

        criteria = Criteria(status='running')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.predicate(self.target))

    def test_status_is_running_with_target_that_is_not_running(self):
        self.backend.status.return_value = Status.UNKNOWN

        criteria = Criteria(status='running')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertFalse(f.predicate(self.target))

    def test_status_is_submitted_with_target_that_is_submitted(self):
        self.backend.status.return_value = Status.SUBMITTED

        criteria = Criteria(status='submitted')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.predicate(self.target))

    def test_status_is_submitted_with_target_that_is_not_submitted(self):
        self.backend.status.return_value = Status.UNKNOWN

        criteria = Criteria(status='submitted')
        f = StatusFilter(self.graph, self.backend, criteria)
        self.assertFalse(f.predicate(self.target))


class TestNameFilter(GWFTestCase):

    def setUp(self):
        self.target1 = Target('TestTarget1', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        self.target2 = Target('TestTarget2', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        self.backend = Mock(name='Backend', spec_set=['status'])
        self.graph = Mock(name='Graph', spec_set=['targets', 'endpoints', 'should_run'])

    def test_should_not_use_filter(self):
        criteria = Criteria(targets=[])
        f = NameFilter(self.graph, self.backend, criteria)
        self.assertFalse(f.use())

    def test_should_use_filter(self):
        criteria = Criteria(targets=['Target1'])
        f = NameFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.use())

    def test_one_name(self):
        self.graph.targets.values.return_value = [self.target1, self.target2]
        criteria = Criteria(targets=['TestTarget1'])
        f = NameFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.predicate(self.target1))
        self.assertFalse(f.predicate(self.target2))

    def test_two_names(self):
        self.graph.targets.values.return_value = [self.target1, self.target2]
        criteria = Criteria(targets=['TestTarget2', 'TestTarget1'])
        f = NameFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.predicate(self.target1))
        self.assertTrue(f.predicate(self.target2))


class TestEndpointFilter(GWFTestCase):

    def setUp(self):
        self.target1 = Target('TestTarget1', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        self.target2 = Target('TestTarget2', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        self.backend = Mock(name='Backend', spec_set=['status'])
        self.graph = Mock(name='Graph', spec_set=['targets', 'endpoints', 'should_run'])

    def test_should_use_filter(self):
        criteria = Criteria(all=False, targets=[])
        f = EndpointFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.use())

    def test_should_not_use_filter(self):
        criteria = Criteria(all=True, targets=[])
        f = EndpointFilter(self.graph, self.backend, criteria)
        self.assertFalse(f.use())

        criteria = Criteria(all=False, targets=['TestTarget1'])
        f = EndpointFilter(self.graph, self.backend, criteria)
        self.assertFalse(f.use())

    def test_with_endpoint_and_non_endpoint_targets(self):
        self.graph.endpoints.return_value = [self.target1]

        criteria = Criteria(all=False, targets=[])
        f = EndpointFilter(self.graph, self.backend, criteria)
        self.assertTrue(f.predicate(self.target1))
        self.assertFalse(f.predicate(self.target2))
