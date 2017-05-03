import logging

from .base import Backend, Status
from ..exceptions import NoLogFoundError

logger = logging.getLogger(__name__)


class TestingBackend(Backend):

    supported_options = set(['cores', 'memory'])
    option_defaults = {'cores': 2, 'memory': '18g'}

    def submit(self, target, dependencies):
        pass

    def cancel(self, target):
        pass

    def status(self, target):
        return Status.UNKNOWN

    def logs(self, target, stderr=False):
        raise NoLogFoundError

    def close(self):
        pass
