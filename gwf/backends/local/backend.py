import json
import logging
import os.path

from .. import Backend
from .client import Client
from .server import State

logger = logging.getLogger(__name__)


class LocalBackend(Backend):
    """A backend that runs targets on a local cluster."""

    supported_options = []
    option_defaults = {}

    def get_client(self):
        return Client(('localhost', int(self.config['workers_port'])))

    def setup_argument_parser(self, parser, subparsers):
        parser.add_argument(
            '--workers-port',
            default=12345,
            type=int,
            help='Port where worker manager accepts clients.',
        )

    def configure(self, *args, **kwargs):
        super().configure(*args, **kwargs)

        self._db_path = os.path.join(
            self.workflow.working_dir,
            '.gwf/local-backend-jobdb.json'
        )

        try:
            with open(self._db_path) as fileobj:
                self._job_db = json.load(fileobj)
        except OSError:
            self._job_db = {}

    def submit(self, target):
        deps = [
            self._job_db[dep.name]
            for dep in self.workflow.dependencies[target]
        ]
        job_id = self.get_client().submit(target, deps=deps)
        self._job_db[target.name] = job_id

    def cancel(self, target):
        pass

    def submitted(self, target):
        return target.name in self._job_db

    def running(self, target):
        job_id = self._job_db[target.name]
        return self.get_client().status()[job_id] == State.running

    def logs(self, target, stderr=False):
        pass

    def close(self):
        try:
            with open(self._db_path, 'w') as fileobj:
                json.dump(self._job_db, fileobj)
        except:
            logger.error(
                'Closing backend caused an exception to be raised.',
                exc_info=True,
            )
