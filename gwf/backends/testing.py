import logging

from . import Backend, Status
from .logmanager import MemoryLogManager

logger = logging.getLogger(__name__)


class TestingBackend(Backend):

    log_manager = MemoryLogManager()

    option_defaults = {'cores': 2, 'memory': '18g', 'nodes': None}

    def submit(self, target, dependencies):
        pass

    def cancel(self, target):
        pass

    def status(self, target):
        return Status.UNKNOWN

    def close(self):
        pass
