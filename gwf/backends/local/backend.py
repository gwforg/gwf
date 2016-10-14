import json
import logging
import os.path

from .. import Backend, FileLogsMixin
from .client import Client
from .server import State

logger = logging.getLogger(__name__)


class LocalBackend(FileLogsMixin, Backend):
    """Backend that runs targets on a local cluster.

    To use this backend you must activate the `local` backend and start a
    local cluster (with one or more workers) that the backend can submit targets
    to. To start a cluster with two workers run the command::

        gwf -b local workers -n 2

    If the local backend is your default backend you can of course omit the
    ``-b local`` option.
    """

    supported_options = []
    option_defaults = {}

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
            self.working_dir,
            '.gwf/local-backend-jobdb.json'
        )

        try:
            with open(self._db_path) as fileobj:
                self._job_db = json.load(fileobj)
        except OSError:
            self._job_db = {}

        self.client = Client(('localhost', self.config['workers_port']))

        status = self.client.status()
        self._job_db = {
            k: v
            for k, v in self._job_db.items()
            if v in status and status[v] != State.completed
        }

    def submit(self, target, dependencies):
        deps = [
            self._job_db[dep.name]
            for dep in dependencies
            if dep.name in self._job_db
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

    def failed(self, target):
        job_id = self._job_db[target.name]
        return self.client.status()[job_id] == State.failed

    def completed(self, target):
        job_id = self._job_db[target.name]
        return self.client.status()[job_id] == State.completed

    def close(self):
        try:
            with open(self._db_path, 'w') as fileobj:
                json.dump(self._job_db, fileobj)
        except:
            logger.error(
                'Closing backend caused an exception to be raised.',
                exc_info=True,
            )
