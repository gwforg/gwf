from gwf.backends import Backend, Status
from gwf.exceptions import LogNotFoundError


class MockBackend(Backend):

    def __init__(self):
        super().__init__()
        self._tracked = {}

    def submit(self, target, dependencies):
        self._tracked[target] = Status.SUBMITTED

    def cancel(self, target):
        del self._tracked[target]

    def status(self, target):
        return self._tracked.get(target, Status.UNKNOWN)

    def logs(self, target, stderr=False):
        raise LogNotFoundError

    def close(self):
        pass

    def set_status(self, target, status):
        assert status in (Status.RUNNING, Status.UNKNOWN)
        self._tracked[target] = status
