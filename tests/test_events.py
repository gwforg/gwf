import unittest
from unittest.mock import Mock

from gwf.events import Event


class TestEvent(unittest.TestCase):

    def setUp(self):
        self.event = Event(name='test_event')
        self.mock_callback = Mock(name='test_callback')

    def test_registering_callback_and_triggering_event_calls_callback(self):
        self.event.register(self.mock_callback)
        self.event.trigger()
        self.mock_callback.assert_called_once_with()

    def test_registering_callback_and_triggering_event_with_args_calls_callback_with_args(self):
        self.event.register(self.mock_callback)
        self.event.trigger('hello', who='world')
        self.mock_callback.assert_called_once_with('hello', who='world')

    def test_registering_callback_twice_and_triggering_does_not_run_it_twice(self):
        self.event.register(self.mock_callback)
        self.event.register(self.mock_callback)
        self.event.trigger()
        self.mock_callback.assert_called_once_with()

    def test_callback_is_not_called_when_it_has_been_unregistered(self):
        self.event.register(self.mock_callback)
        self.event.trigger()
        self.event.unregister(self.mock_callback)
        self.event.trigger()
        self.mock_callback.assert_called_once_with()
