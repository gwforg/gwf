"""Provides easy, dynamic filtering of targets."""

from .backends.base import Status


FILTERS = []


def register_filter(filter_cls):
    FILTERS.append(filter_cls)


class Criteria:
    """A container for filtering criteria."""
    def __init__(self, **kwargs):
        self.__dict__ = kwargs


class FilterType(type):
    """A metaclass for filters.

    The only purpose of this metaclass is to automatically register filters.
    """

    def __new__(mcs, name, bases, class_dict):
        cls = type.__new__(mcs, name, bases, class_dict)
        if bases:
            register_filter(cls)
        return cls


class Filter(metaclass=FilterType):
    """A base class for filters."""
    def __init__(self, graph, backend, criteria):
        self.graph = graph
        self.backend = backend
        self.criteria = criteria

    def apply(self, targets):
        return (target for target in targets if self.predicate(target))


class StatusFilter(Filter):

    def use(self):
        return self.criteria.status

    def predicate(self, target):
        if self.criteria.status == 'completed':
            return not self.graph.should_run(target)
        if self.criteria.status == 'shouldrun':
            return self.graph.should_run(target) and self.backend.status(target) == Status.UNKNOWN
        if self.criteria.status == 'running':
            return self.backend.status(target) == Status.RUNNING
        if self.criteria.status == 'submitted':
            return self.backend.status(target) == Status.SUBMITTED


class NameFilter(Filter):
    def use(self):
        return self.criteria.targets

    def predicate(self, target):
        return target.name in self.criteria.targets


class EndpointFilter(Filter):
    def use(self):
        return not self.criteria.all and not self.criteria.targets

    def predicate(self, target):
        return target in self.graph.endpoints()


def filter(graph, backend, criteria):
    targets = iter(graph.targets.values())
    for filter_cls in FILTERS:
        filter = filter_cls(graph, backend, criteria)
        if filter.use():
            targets = filter.apply(targets)
    return targets
