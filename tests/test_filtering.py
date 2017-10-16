from unittest.mock import create_autospec, Mock

import pytest

from gwf import Target
from gwf.core import Graph, Scheduler
from gwf.backends import Status
from gwf.backends.testing import TestingBackend
from gwf.filtering import StatusFilter, Criteria, NameFilter, EndpointFilter, filter


@pytest.fixture
def backend():
    return create_autospec(TestingBackend(), spec_set=True)


@pytest.fixture
def graph():
    return create_autospec(Graph(dependencies={}, dependents={}, provides={}, targets={}), spec_set=True)


@pytest.fixture
def scheduler(backend, graph):
    scheduler = Scheduler(backend=backend, graph=graph)
    scheduler.should_run = create_autospec(scheduler.should_run, spec_set=True)
    return scheduler


@pytest.fixture
def status_filter(scheduler):
    return StatusFilter(scheduler=scheduler)


@pytest.fixture
def name_filter():
    return NameFilter()


@pytest.fixture
def endpoint_filter():
    return EndpointFilter(endpoints=set())


def test_filter_status(status_filter):
    criteria = Criteria(status='completed')
    assert status_filter.use(criteria)

    criteria = Criteria(status=None)
    assert not status_filter.use(criteria)


def test_filter_status_completed(status_filter):
    criteria = Criteria(status='completed')
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

    status_filter.scheduler.should_run.return_value = False
    assert list(status_filter.apply([target], criteria)) == [target]

    status_filter.scheduler.should_run.return_value = True
    assert list(status_filter.apply([target], criteria)) == []


def test_filter_status_shouldrun(status_filter):
    criteria = Criteria(status='shouldrun')
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

    status_filter.scheduler.backend.status.return_value = Status.UNKNOWN
    status_filter.scheduler.should_run.return_value = False
    assert list(status_filter.apply([target], criteria)) == []

    status_filter.scheduler.backend.status.return_value = Status.UNKNOWN
    status_filter.scheduler.should_run.return_value = True
    assert list(status_filter.apply([target], criteria)) == [target]

    status_filter.scheduler.backend.status.return_value = Status.RUNNING
    status_filter.scheduler.should_run.return_value = False
    assert list(status_filter.apply([target], criteria)) == []

    status_filter.scheduler.backend.status.return_value = Status.RUNNING
    status_filter.scheduler.should_run.return_value = True
    assert list(status_filter.apply([target], criteria)) == []

    status_filter.scheduler.backend.status.return_value = Status.SUBMITTED
    status_filter.scheduler.should_run.return_value = False
    assert list(status_filter.apply([target], criteria)) == []

    status_filter.scheduler.backend.status.return_value = Status.SUBMITTED
    status_filter.scheduler.should_run.return_value = True
    assert list(status_filter.apply([target], criteria)) == []


def test_filter_status_running(status_filter):
    criteria = Criteria(status='running')
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

    status_filter.scheduler.backend.status.return_value = Status.UNKNOWN
    assert list(status_filter.apply([target], criteria)) == []

    status_filter.scheduler.backend.status.return_value = Status.SUBMITTED
    assert list(status_filter.apply([target], criteria)) == []

    status_filter.scheduler.backend.status.return_value = Status.RUNNING
    assert list(status_filter.apply([target], criteria)) == [target]


def test_filter_status_submitted(status_filter):
    criteria = Criteria(status='submitted')
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

    status_filter.scheduler.backend.status.return_value = Status.UNKNOWN
    assert list(status_filter.apply([target], criteria)) == []

    status_filter.scheduler.backend.status.return_value = Status.SUBMITTED
    assert list(status_filter.apply([target], criteria)) == [target]

    status_filter.scheduler.backend.status.return_value = Status.RUNNING
    assert list(status_filter.apply([target], criteria)) == []


def test_filter_name(name_filter):
    criteria = Criteria(targets=['TestTarget'])
    assert name_filter.use(criteria)

    criteria = Criteria(targets=[])
    assert not name_filter.use(criteria)

    target1 = Target('Foo', inputs=[], outputs=[], options={}, working_dir='/some/dir')
    target2 = Target('Bar', inputs=[], outputs=[], options={}, working_dir='/some/dir')
    target3 = Target('FooBar', inputs=[], outputs=[], options={}, working_dir='/some/dir')

    criteria = Criteria(targets=['Foo'])
    assert set(name_filter.apply([target1, target2, target3], criteria)) == {target1}

    criteria = Criteria(targets=['Foo*'])
    assert set(name_filter.apply([target1, target2, target3], criteria)) == {target1, target3}

    criteria = Criteria(targets=['Foo', 'Bar'])
    assert set(name_filter.apply([target1, target2, target3], criteria)) == {target1, target2}


def test_filter_endpoint(endpoint_filter):
    criteria = Criteria(all=True, targets=['Foo'])
    assert not endpoint_filter.use(criteria)

    criteria = Criteria(all=False, targets=['Foo'])
    assert not endpoint_filter.use(criteria)

    criteria = Criteria(all=True, targets=[])
    assert not endpoint_filter.use(criteria)

    criteria = Criteria(all=False, targets=[])
    assert endpoint_filter.use(criteria)

    target1 = Target('Foo', inputs=[], outputs=[], options={}, working_dir='/some/dir')
    target2 = Target('Bar', inputs=[], outputs=[], options={}, working_dir='/some/dir')
    target3 = Target('FooBar', inputs=[], outputs=[], options={}, working_dir='/some/dir')

    endpoint_filter.endpoints = {target1}
    assert set(endpoint_filter.apply([target1, target2, target3], criteria)) == {target1}

    endpoint_filter.endpoints = {target1, target3}
    assert set(endpoint_filter.apply([target1, target2, target3], criteria)) == {target1, target3}
