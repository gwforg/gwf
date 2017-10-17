import fnmatch

from .backends import Status


class Criteria:
    """A container for filtering criteria."""
    def __init__(self, **kwargs):
        self.__dict__ = kwargs


class Filter:
    """A base class for filters."""

    def apply(self, targets, criteria):
        return (target for target in targets if self.predicate(target, criteria))

    def predicate(self, target, criteria):
        return True


class StatusFilter(Filter):

    def __init__(self, scheduler):
        self.scheduler = scheduler

    def use(self, criteria):
        return criteria.status

    def predicate(self, target, criteria):
        if criteria.status == 'completed':
            return not self.scheduler.should_run(target)
        if criteria.status == 'shouldrun':
            return (
                self.scheduler.should_run(target) and
                self.scheduler.backend.status(target) == Status.UNKNOWN
            )
        if criteria.status == 'running':
            return self.scheduler.backend.status(target) == Status.RUNNING
        if criteria.status == 'submitted':
            return self.scheduler.backend.status(target) == Status.SUBMITTED


class NameFilter(Filter):

    def use(self, criteria):
        return criteria.targets

    def apply(self, targets, criteria):
        target_name_map = {target.name: target for target in targets}
        return {
            target_name_map[name]
            for pattern in criteria.targets
            for name in fnmatch.filter(target_name_map.keys(), pattern)
        }


class EndpointFilter(Filter):

    def __init__(self, endpoints, negate=False):
        self.endpoints = endpoints
        self.negate = negate

    def use(self, criteria):
        return not criteria.all and not criteria.targets

    def predicate(self, target, criteria):
        if self.negate:
            return target not in self.endpoints
        return target in self.endpoints


class CompositeFilter(Filter):

    def __init__(self, filters=None):
        self.filters = filters

    def apply(self, targets, criteria):
        for filter in self.filters:
            if filter.use(criteria):
                targets = filter.apply(targets, criteria)
        return targets


def filter_generic(targets, criteria, filters):
    """Filter targets given some criteria`and filters.

    Return all targets from `targets` matching the criteria in `criteria` using `filters`.
    """
    filterer = CompositeFilter(filters)
    return filterer.apply(targets, criteria)


def filter_names(targets, patterns):
    """Filter targets with a list of patterns.

    Return all targets in `targets` where the target name matches one or more of the patterns in `pattern`.
    """
    filter = NameFilter()
    criteria = Criteria(targets=patterns)
    if filter.use(criteria):
        return filter.apply(targets, criteria=criteria)
    return targets