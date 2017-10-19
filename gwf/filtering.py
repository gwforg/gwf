import fnmatch

from .backends import Status


class ApplyMixin:
    """A mixin for predicate-based filters providing the `apply` method.

    Most filters are predicate-based in the sense that they simply filter targets one by one based on a predicate
    function that decides whether to include the target or not. Such filters can inherit this mixin and then only need
    to declare a :func:`predicate` method which returns `True` if the target should be included and `False` otherwise.

    For examples of using this mixin, see the :class:`StatusFilter` and :class:`EndpointFilter` filters.
    """

    def apply(self, targets):
        """Apply the filter to all `targets`.

        This method returns a generator yielding all targets in `targets` for each :func:`predicate` returns `True`.
        """
        return (target for target in targets if self.predicate(target))

    def predicate(self, target):
        """Return `True` if `target` should be included, `False` otherwise.

        This method must be overriden by subclasses.
        """
        raise NotImplementedError('predicate() must be implemented by subclasses.')


class StatusFilter(ApplyMixin):

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


class NameFilter:

    def __init__(self, patterns):
        self.patterns = patterns

    def apply(self, targets):
        target_name_map = {target.name: target for target in targets}
        return {
            target_name_map[name]
            for pattern in self.patterns
            for name in fnmatch.filter(target_name_map.keys(), pattern)
        }


class EndpointFilter(ApplyMixin):

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


class CompositeFilter:

    def __init__(self, filters=None):
        self.filters = filters

    def apply(self, targets):
        for filter in self.filters:
            targets = filter.apply(targets)
        return targets


def filter_generic(targets, filters):
    """Filter targets given a list of filters.

    Return all targets from `targets` passing all `filters`. For example:

    .. code-block:: python

        matched_targets = filter_generic(
            targets=graph.targets.values(),
            filters=[
                NameFilter(patterns=['Foo*'],
                StatusFilter(scheduler=scheduler, status='running'),
            ]
        )

    returns a generator yielding all targets with a name matching ``Foo*`` which are currently running.

    :arg targets:
        A list of targets to be filtered.

    :arg filters:
        A list of :class:`Filter` instances.
    """
    return CompositeFilter(filters).apply(targets)


def filter_names(targets, patterns):
    """Filter targets with a list of patterns.

    Return all targets in `targets` where the target name matches one or more of the patterns in `pattern`. For example:

    .. code-block:: python

        matched_targets = filter_names(graph.targets.values(), ['Foo*'])

    returns a generator yielding all targets with a name matching the pattern `Foo*`. Multiple patterns can be provided:

    .. code-block:: python

        matched_targets = filter_names(graph.targets.values(), ['Foo*', 'Bar*'])

    returns all targets with a name matching either ``Foo*`` or ``Bar*``.

    This function is a simple wrapper around :class:`NameFilter`.
    """
    return NameFilter(patterns=patterns).apply(targets)
