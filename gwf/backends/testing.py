import logging

from . import Backend

logger = logging.getLogger(__name__)


class TestingBackend(Backend):

    name = 'testing'
    supported_options = set(['cores', 'memory'])

    def configure(self, *args, **kwargs):
        pass

    def submit(self, target):
        pass

    def cancel(self, target):
        pass

    def submitted(self, target):
        return False

    def running(self, target):
        return False
