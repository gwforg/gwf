import logging

logger = logging.getLogger(__name__)


class Event:
    """A simple event implementation.

    This class implements the observer pattern in a reusable fashion. Event
    objects can be created anywhere (in modules, class definitions etc.) and
    anyone can register for the event or trigger it.

    .. note::
        Even though it is possible to trigger any event, please only trigger
        your own events unless you really know what you're doing.

    Triggering an event will call all registered callbacks with the arguments
    given to :func:`trigger`.
    """

    def __init__(self, name):
        self.name = name
        self._callbacks = set()
        self._logger = logging.getLogger('{}.{}'.format(__name__, self.name))

    def register(self, callback):
        """Register a callback to this event.

        :param callable callback:
            a callable do call whenever the event is triggered.

        Registering the same callback twice will only add the callback once.
        """
        if callback not in self._callbacks:
            self._logger.debug('Registered callback: %r.', callback)
            self._callbacks.add(callback)

    def unregister(self, callback):
        """Unregister a callback.

        :param callable callback:
            a callable to unregister.

        It is safe to unregister the same callback twice. The second call to
        :func:`unregister` will be a no-op.
        """
        self._logger.debug('Unregistered callback: %r.', callback)
        self._callbacks.discard(callback)

    def trigger(self, *args, **kwargs):
        """Trigger the event.

        All arguments will be forwarded to the callbacks. Callbacks are called
        in no particular order.
        """
        for callback in self._callbacks:
            self._logger.debug(
                'Triggering callback: %r with args %r and kwargs %r.',
                callback,
                args,
                kwargs,
            )
            callback(*args, **kwargs)

    def __repr__(self):
        return '{}(name={})'.format(self.__class__.__name__, self.name)


pre_schedule = Event('pre_schedule')
post_schedule = Event('post_schedule')
