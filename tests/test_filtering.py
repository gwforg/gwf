import pytest

from gwf import Target
from gwf.core import TargetStatus
from gwf.filtering import StatusFilter, NameFilter, EndpointFilter


def make_status_provider(statuses):
    def status(target):
        return statuses[target]

    return status


target1 = Target(
    "TestTarget1", inputs=[], outputs=[], options={}, working_dir="/some/dir"
)
target2 = Target(
    "TestTarget2", inputs=[], outputs=[], options={}, working_dir="/some/dir"
)
target3 = Target(
    "TestTarget3", inputs=[], outputs=[], options={}, working_dir="/some/dir"
)
target4 = Target(
    "TestTarget4", inputs=[], outputs=[], options={}, working_dir="/some/dir"
)
target5 = Target(
    "TestTarget5", inputs=[], outputs=[], options={}, working_dir="/some/dir"
)
target6 = Target(
    "TestTarget6", inputs=[], outputs=[], options={}, working_dir="/some/dir"
)

targets = [target1, target2, target3, target4, target5, target6]

status_provider = make_status_provider(
    {
        target1: TargetStatus.COMPLETED,
        target2: TargetStatus.RUNNING,
        target3: TargetStatus.SHOULDRUN,
        target4: TargetStatus.SUBMITTED,
        target5: TargetStatus.SHOULDRUN,
        target6: TargetStatus.COMPLETED,
    }
)


@pytest.mark.parametrize(
    "status,result",
    [
        ([TargetStatus.COMPLETED], [target1, target6]),
        ([TargetStatus.RUNNING], [target2]),
        ([TargetStatus.SHOULDRUN], [target3, target5]),
        ([TargetStatus.SUBMITTED], [target4]),
        ([TargetStatus.SUBMITTED, TargetStatus.RUNNING], [target2, target4]),
        (
            [TargetStatus.COMPLETED, TargetStatus.SHOULDRUN],
            [target1, target6, target3, target5],
        ),
    ],
)
def test_filter_status(status, result):
    status_filter = StatusFilter(status_provider=status_provider, status=status,)
    assert set(status_filter.apply(targets)) == set(result)


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
