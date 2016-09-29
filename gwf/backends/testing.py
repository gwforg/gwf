import logging

from . import Backend

logger = logging.getLogger(__name__)


class TestingBackend(Backend):

    name = 'testing'
    supported_options = set(['cores', 'memory'])
    option_defaults = {'cores': 2, 'memory': '18gb'}

    def submit(self, target):
        pass

    def cancel(self, target):
        pass

    def submitted(self, target):
        return False

    def running(self, target):
        return False

    def logs(self, target, stderr=False, rewind=0):
        pass
