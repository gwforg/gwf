from unittest.mock import create_autospec

import pytest

from gwf import Target
from gwf.core import Graph, Scheduler, TargetStatus
from gwf.backends import Status
from gwf.backends.testing import TestingBackend
from gwf.filtering import StatusFilter, NameFilter, EndpointFilter


@pytest.fixture
def backend():
    return create_autospec(TestingBackend(), spec_set=True)


@pytest.fixture
def graph():
    return create_autospec(
        Graph(dependencies={}, dependents={}, provides={}, targets={}, unresolved={}),
        spec_set=True,
    )


@pytest.fixture
def scheduler(backend, graph):
    scheduler = Scheduler(backend=backend, graph=graph)
    scheduler.should_run = create_autospec(scheduler.should_run, spec_set=True)
    return scheduler


def test_filter_status_completed(scheduler):
    scheduler.backend.status.return_value = Status.UNKNOWN

    status_filter = StatusFilter(scheduler=scheduler, status=[TargetStatus.COMPLETED])
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    status_filter.scheduler.should_run.return_value = False
    assert list(status_filter.apply([target])) == [target]

    status_filter.scheduler.should_run.return_value = True
    assert list(status_filter.apply([target])) == []


def test_filter_status_shouldrun(scheduler):
    status_filter = StatusFilter(scheduler=scheduler, status=[TargetStatus.SHOULDRUN])
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    status_filter.scheduler.backend.status.return_value = Status.UNKNOWN
    status_filter.scheduler.should_run.return_value = False
    assert list(status_filter.apply([target])) == []

    status_filter.scheduler.backend.status.return_value = Status.UNKNOWN
    status_filter.scheduler.should_run.return_value = True
    assert list(status_filter.apply([target])) == [target]

    status_filter.scheduler.backend.status.return_value = Status.RUNNING
    status_filter.scheduler.should_run.return_value = False
    assert list(status_filter.apply([target])) == []

    status_filter.scheduler.backend.status.return_value = Status.RUNNING
    status_filter.scheduler.should_run.return_value = True
    assert list(status_filter.apply([target])) == []

    status_filter.scheduler.backend.status.return_value = Status.SUBMITTED
    status_filter.scheduler.should_run.return_value = False
    assert list(status_filter.apply([target])) == []

    status_filter.scheduler.backend.status.return_value = Status.SUBMITTED
    status_filter.scheduler.should_run.return_value = True
    assert list(status_filter.apply([target])) == []


def test_filter_status_running(scheduler):
    status_filter = StatusFilter(scheduler=scheduler, status=[TargetStatus.RUNNING])
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    status_filter.scheduler.backend.status.return_value = Status.UNKNOWN
    assert list(status_filter.apply([target])) == []

    status_filter.scheduler.backend.status.return_value = Status.SUBMITTED
    assert list(status_filter.apply([target])) == []

    status_filter.scheduler.backend.status.return_value = Status.RUNNING
    assert list(status_filter.apply([target])) == [target]


def test_filter_status_submitted(scheduler):
    status_filter = StatusFilter(scheduler=scheduler, status=[TargetStatus.SUBMITTED])
    target = Target(
        "TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    status_filter.scheduler.backend.status.return_value = Status.UNKNOWN
    assert list(status_filter.apply([target])) == []

    status_filter.scheduler.backend.status.return_value = Status.SUBMITTED
    assert list(status_filter.apply([target])) == [target]

    status_filter.scheduler.backend.status.return_value = Status.RUNNING
    assert list(status_filter.apply([target])) == []


def test_filter_name():
    target1 = Target("Foo", inputs=[], outputs=[], options={}, working_dir="/some/dir")
    target2 = Target("Bar", inputs=[], outputs=[], options={}, working_dir="/some/dir")
    target3 = Target(
        "FooBar", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    name_filter = NameFilter(patterns=["Foo"])
    assert set(name_filter.apply([target1, target2, target3])) == {target1}

    name_filter = NameFilter(patterns=["Foo*"])
    assert set(name_filter.apply([target1, target2, target3])) == {target1, target3}

    name_filter = NameFilter(patterns=["Foo", "Bar"])
    assert set(name_filter.apply([target1, target2, target3])) == {target1, target2}


def test_filter_endpoint():
    target1 = Target("Foo", inputs=[], outputs=[], options={}, working_dir="/some/dir")
    target2 = Target("Bar", inputs=[], outputs=[], options={}, working_dir="/some/dir")
    target3 = Target(
        "FooBar", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )

    endpoint_filter = EndpointFilter(endpoints={target1})
    assert set(endpoint_filter.apply([target1, target2, target3])) == {target1}

    endpoint_filter = EndpointFilter(endpoints={target1, target3})
    assert set(endpoint_filter.apply([target1, target2, target3])) == {target1, target3}

    endpoint_filter = EndpointFilter(endpoints={target1}, mode="exclude")
    assert set(endpoint_filter.apply([target1, target2, target3])) == {target2, target3}

    endpoint_filter = EndpointFilter(endpoints={target1, target3}, mode="exclude")
    assert set(endpoint_filter.apply([target1, target2, target3])) == {target2}
