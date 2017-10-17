import fnmatch

from .backends import Status


class Criteria:
    """A container for filtering criteria."""
    def __init__(self, **kwargs):
        self.__dict__ = kwargs


class Filter:
    """A base class for filters."""

    def apply(self, targets):
        return (target for target in targets if self.predicate(target))

    def predicate(self, target):
        return True


class StatusFilter(Filter):

    def __init__(self, scheduler, status):
        self.scheduler = scheduler
        self.status = status

    def predicate(self, target):
        if self.status == 'completed':
            return not self.scheduler.should_run(target)
        if self.status == 'shouldrun':
            return (
                self.scheduler.should_run(target) and
                self.scheduler.backend.status(target) == Status.UNKNOWN
            )
        if self.status == 'running':
            return self.scheduler.backend.status(target) == Status.RUNNING
        if self.status == 'submitted':
            return self.scheduler.backend.status(target) == Status.SUBMITTED


class NameFilter(Filter):

    def __init__(self, patterns):
        self.patterns = patterns

    def apply(self, targets):
        target_name_map = {target.name: target for target in targets}
        return {
            target_name_map[name]
            for pattern in self.patterns
            for name in fnmatch.filter(target_name_map.keys(), pattern)
        }


class EndpointFilter(Filter):

    def __init__(self, endpoints, mode='include'):
        self.endpoints = endpoints
        self.mode = mode

    def predicate(self, target):
        if self.mode == 'exclude':
            return target not in self.endpoints
        elif self.mode == 'include':
            return target in self.endpoints
        else:
            raise ValueError('Argument mode must be either "include" or "exclude".')


class CompositeFilter(Filter):

    def __init__(self, filters=None):
        self.filters = filters

    def apply(self, targets):
        for filter in self.filters:
            targets = filter.apply(targets)
        return targets


def filter_generic(targets, filters):
    """Filter targets given a list of filters.

    :arg targets:
        A list of targets to be filtered. This will most often be all of the targets in a graph.

    :arg filters:
        A list of :class:`Filter` instances.
    """
    return CompositeFilter(filters).apply(targets)


def filter_names(targets, patterns):
    """Filter targets with a list of patterns.

    Return all targets in `targets` where the target name matches one or more of the patterns in `pattern`.
    """
    return NameFilter(patterns=patterns).apply(targets)
