import logging
import multiprocessing
import os
import os.path
import subprocess
import sys
import uuid
from enum import Enum
from multiprocessing import Manager
from multiprocessing.connection import Client as Client_
from multiprocessing.connection import Listener
from multiprocessing.pool import Pool

from . import Backend, Status
from .exceptions import BackendError, UnsupportedOperationError, UnknownDependencyError
from .logmanager import FileLogManager
from ..conf import config
from ..utils import PersistableDict

__all__ = ('Client', 'Server', 'LocalBackend',)

logger = logging.getLogger(__name__)


def _gen_task_id():
    return uuid.uuid4().hex


def catch_keyboard_interrupt(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except KeyboardInterrupt:
            logger.debug('Shutting down worker %s...', multiprocessing.current_process().pid)
            sys.exit(0)

    return wrapper


class ServerError(BackendError):
    pass


class LocalStatus(Enum):
    UNKNOWN = 0
    SUBMITTED = 1
    RUNNING = 2
    FAILED = 3
    COMPLETED = 4


class Request:
    def handle(self, task_queue, status_queue):
        """Handle this request."""


class SubmitRequest(Request):
    def __init__(self, target, deps, stdout_path, stderr_path):
        self.target = target
        self.deps = deps
        self.stdout_path = stdout_path
        self.stderr_path = stderr_path

    def handle(self, task_queue, status_dict):
        task_id = _gen_task_id()
        status_dict[task_id] = LocalStatus.SUBMITTED
        task_queue.put((task_id, self))
        logger.debug('Task %s was queued with id %s.', self.target.name, task_id)
        return task_id


class StatusRequest(Request):
    def handle(self, task_queue, status_dict):
        return dict(status_dict)


class Client:
    """A client for communicating with the workers."""

    def __init__(self, *args, **kwargs):
        self.client = Client_(*args, **kwargs)

    def submit(self, target, stdout_path, stderr_path, deps=None):
        if deps is None:
            deps = []
        request = SubmitRequest(target=target, deps=deps, stdout_path=stdout_path, stderr_path=stderr_path)
        self.client.send(request)
        return self.client.recv()

    def status(self):
        request = StatusRequest()
        self.client.send(request)
        return self.client.recv()

    def close(self):
        self.client.close()


class LocalBackend(Backend):
    """Backend that runs targets on a local cluster.

    To use this backend you must activate the `local` backend and start a
    local cluster (with one or more workers) that the backend can submit targets
    to. To start a cluster with two workers run the command::

        gwf -b local workers -n 2

    in the working directory of your project. The workflow file must be accessible
    to *gwf*. Thus, if your workflow file is not called `workflow.py` or the
    workflow object is not called `gwf`, you must specify this so that *gwf* can
    locate the workflow::

        gwf -f myworkflow.py:wf1 -b local workers -n 2

    If the local backend is your default backend you can of course omit the
    ``-b local`` option.

    If the ``-n`` option is omitted, *gwf* will detect the number of cores
    available and use all of them.

    To run your workflow, open another terminal and then type::

        gwf -b local run

    To stop the pool of workers press :kbd:`Control-c`.

    **Backend options:**

    * **local.host (str):** Set the host that the workers are running on (default: localhost).
    * **local.port (int):** Set the port used to connect to the workers (default: 12345).

    **Target options:**

    None available.
    """

    log_manager = FileLogManager()

    option_defaults = {}

    def __init__(self):
        super().__init__()

        self._tracked = PersistableDict(os.path.join('.gwf/local-backend-tracked.json'))

        host = config.get('local.host', 'localhost')
        port = config.get('local.port', 12345)
        try:
            self.client = Client((host, port))
        except ConnectionRefusedError:
            raise BackendError(
                'Local backend could not connect to workers on port {}. '
                'Workers can be started by running "gwf workers". '
                'You can read more in the documentation: '
                'http://gwf.readthedocs.io/en/latest/reference/backends.html#local'.format(port)
            )

        self._status = self.client.status()
        for target_name, target_job_id in list(self._tracked.items()):
            if target_job_id not in self._status or self._status[target_job_id] == LocalStatus.COMPLETED:
                del self._tracked[target_name]

    def submit(self, target, dependencies):
        try:
            dependency_ids = [self._tracked[dep.name] for dep in dependencies]
        except KeyError as exc:
            key, = exc.args
            raise UnknownDependencyError(key)

        task_id = self.client.submit(
            target,
            deps=dependency_ids,
            stdout_path=self.log_manager.stdout_path(target),
            stderr_path=self.log_manager.stderr_path(target)
        )
        self._tracked[target.name] = task_id
        self._status[task_id] = LocalStatus.SUBMITTED

    def cancel(self, target):
        raise UnsupportedOperationError('cancel')

    def status(self, target):
        try:
            target_job_id = self._tracked[target.name]
            target_status = self._status[target_job_id]
            if target_status == LocalStatus.RUNNING:
                return Status.RUNNING
            elif target_status == LocalStatus.SUBMITTED:
                return Status.SUBMITTED
            else:
                return Status.UNKNOWN
        except KeyError:
            return Status.UNKNOWN

    def close(self):
        self._tracked.persist()


class Worker:
    def __init__(self, status, queue, waiting):
        self.status = status
        self.queue = queue
        self.waiting = waiting

        self.run()

    @catch_keyboard_interrupt
    def run(self):
        while True:
            task_id, request = self.queue.get()

            # If the task isn't pending, it may have been resubmitted by multiple
            # tasks in different workers. We shouldn't run it twice, so we'll skip
            # it.
            if self.status[task_id] != LocalStatus.SUBMITTED:
                continue

            if not self.check_dependencies(task_id, request):
                continue

            self.handle_task(task_id, request)

    def handle_task(self, task_id, request):
        logger.debug('Task %s started target %r.', task_id, request.target)
        self.status[task_id] = LocalStatus.RUNNING

        try:
            self.execute_target(request.target, stdout_path=request.stdout_path, stderr_path=request.stderr_path)
        except:
            self.status[task_id] = LocalStatus.FAILED
            logger.error('Task %s failed.', task_id, exc_info=True)
        else:
            self.status[task_id] = LocalStatus.COMPLETED
            logger.debug('Task %s completed target %r.', task_id, request.target)
        finally:
            self.requeue_dependents(task_id)

    def check_dependencies(self, task_id, request):
        """Check dependencies before running a task.

        :return: `True` if the task can run, `False` if not."""
        any_dep_failed = any(
            self.status[dep_id] == LocalStatus.FAILED
            for dep_id in request.deps
        )

        if any_dep_failed:
            self.status[task_id] = LocalStatus.FAILED
            logger.error(
                'Task %s failed since a dependency failed.',
                task_id,
                exc_info=True
            )
            return False

        has_non_satisfied_dep = False
        for dep_id in request.deps:
            if self.status[dep_id] != LocalStatus.COMPLETED:
                logger.debug('Task %s set to wait for %s.', task_id, dep_id)

                if dep_id not in self.waiting:
                    self.waiting[dep_id] = []

                # Do this weird thing to sync the dictionary. Using
                # deps_dict[dep_id].append(...) doesn't work.
                tmp = self.waiting[dep_id]
                tmp.append((task_id, request))
                self.waiting[dep_id] = tmp

                has_non_satisfied_dep = True

        if has_non_satisfied_dep:
            return False
        return True

    def requeue_dependents(self, task_id):
        """Requeue targets that waiting on a specific task."""
        if task_id not in self.waiting:
            return

        logger.debug('Task %s has waiting dependents. Requeueing.', task_id)
        for dep_task_id, dep_request in self.waiting[task_id]:
            self.queue.put((dep_task_id, dep_request))

    def execute_target(self, target, stdout_path, stderr_path):
        env = os.environ.copy()
        env['GWF_TARGET_NAME'] = target.name

        with open(stdout_path, mode='w') as stdout_fp, open(stderr_path, mode='w') as stderr_fp:
            process = subprocess.Popen(
                ['bash'],
                stdin=subprocess.PIPE,
                stdout=stdout_fp,
                stderr=stderr_fp,
                universal_newlines=True,
                cwd=target.working_dir,
                env=env,
            )

            process.communicate(target.spec)
            if process.returncode != 0:
                raise Exception('Target {} exited with a non-zero return code.'.format(target.name))


class Server:
    def __init__(self, hostname='', port=0, num_workers=None):
        self.hostname = hostname
        self.port = port
        self.num_workers = num_workers

        self.manager = Manager()
        self.status = self.manager.dict()
        self.queue = self.manager.Queue()
        self.waiting = self.manager.dict()

    def handle_request(self, request):
        try:
            logger.debug('Received request %r.', request)
            return request.handle(self.queue, self.status)
        except:
            logger.error('Invalid request %r.', request, exc_info=True)

    def handle_client(self, conn):
        logger.debug('Accepted client connection.')
        try:
            while True:
                request = conn.recv()
                response = self.handle_request(request)
                if response is not None:
                    conn.send(response)
        except EOFError:
            logger.debug('Client connection closed.')

    def wait_for_clients(self, serv):
        while True:
            client = serv.accept()
            self.handle_client(client)

    def start(self):
        """Starts a server that controls local workers.

        Calling this function starts a pool of `num_workers` workers used to run
        targets sent to the server. The server will run indefinitely unless shut
        down by the user.
        """
        try:
            serv = Listener((self.hostname, self.port))
            workers = Pool(
                processes=self.num_workers,
                initializer=Worker,
                initargs=(self.status, self.queue, self.waiting),
            )

            logging.info('Started %s workers, listening on port %s.', self.num_workers, serv.address[1])
            self.wait_for_clients(serv)
        except OSError as e:
            if e.errno == 48:
                raise ServerError(
                    ('Could not start workers listening on port {}. '
                     'The port may already be in use.').format(self.port)
                )
        except KeyboardInterrupt:
            logging.info('Shutting down...')
            workers.close()
            workers.join()
            self.manager.shutdown()
