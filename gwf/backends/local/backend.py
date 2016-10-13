import json
import logging
import os.path

from .. import Backend
from ...exceptions import BackendError
from .client import Client
from .server import State

logger = logging.getLogger(__name__)


class LocalBackend(Backend):
    """A backend that runs targets on a local cluster."""

    name = 'local'

    supported_options = []
    option_defaults = {}

    def configure(self, *args, **kwargs):
        super().configure(*args, **kwargs)
        self.client = Client(('localhost', 25000))

        self._db_path = os.path.join(
            self.workflow.working_dir, '.gwf/local-backend-jobdb.json'
        )
        try:
            with open(self._db_path) as fileobj:
                self._job_db = json.load(fileobj)
        except OSError:
            self._job_db = {}

        status = self.client.status()
        self._job_db = {
            k: v for k, v in self._job_db.items() if v in status
        }

    def submit(self, target):
        deps = [
            self._job_db[dep.name]
            for dep in self.workflow.dependencies[target]
        ]
        job_id = self.client.submit(target, deps=deps)
        self._job_db[target.name] = job_id

    def cancel(self, target):
        pass

    def submitted(self, target):
        return target.name in self._job_db

    def running(self, target):
        job_id = self._job_db[target.name]
        return self.client.status()[job_id] == State.running

    def logs(self, target, stderr=False):
        pass

    def close(self):
        with open(self._db_path, 'w') as fileobj:
            json.dump(self._job_db, fileobj)
        self.client.close()
