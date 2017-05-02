import logging

from .base import Backend, Status
from ..exceptions import NoLogFoundError

logger = logging.getLogger(__name__)


class TestingBackend(Backend):  # pragma: no cover

    supported_options = set(['cores', 'memory'])
    option_defaults = {'cores': 2, 'memory': '18gb'}

    def __init__(self, working_dir):
        self.working_dir = working_dir

    def submit(self, target, dependencies):
        pass

    def cancel(self, target):
        pass

    def status(self, target):
        return Status.UNKNOWN

    def logs(self, target, stderr=False):
        raise NoLogFoundError
