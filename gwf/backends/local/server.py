import logging
import os.path
import subprocess
import sys
import traceback
import uuid
from multiprocessing import Manager
from multiprocessing.connection import Listener
from multiprocessing.pool import ThreadPool as Pool

from ...exceptions import BackendError

logger = logging.getLogger(__name__)


def _gen_task_id():
    return uuid.uuid4().hex


def _get_log_dir_path(target):
    return os.path.join(target.working_dir, '.gwf', 'logs')


def _get_stdout_path(target):
    return os.path.join(
        _get_log_dir_path(target),
        '{}.stdout'.format(target.name)
    )


def _get_stderr_path(target):
    return os.path.join(
        _get_log_dir_path(target),
        '{}.stderr'.format(target.name)
    )


def _make_log_dir(target):
    log_dir = _get_log_dir_path(target)
    try:
        os.makedirs(log_dir)
    except OSError:
        pass
    return log_dir


class ServerError(BackendError):
    pass


class State:
    pending = 0
    started = 1
    completed = 2
    failed = 3


class Request:

    def handle(self, task_queue, status_queue):
        raise NotImplemented()


class SubmitRequest(Request):

    def __init__(self, target, deps):
        self.target = target
        self.deps = deps

    def handle(self, task_queue, status_dict):
        task_id = _gen_task_id()
        status_dict[task_id] = State.pending
        task_queue.put((task_id, self))
        return task_id


class StatusRequest(Request):

    def handle(self, task_queue, status_dict):
        return dict(status_dict)


def worker(task_queue, status_dict, deps_dict, deps_lock):
    while True:
        task_id, request = task_queue.get()

        # If the task isn't pending, it may have been resubmitted by multiple
        # tasks in different workers. We shouldn't run it twice, so we'll skip
        # it.
        if status_dict[task_id] != State.pending:
            continue

        with deps_lock:
            any_dep_failed = any(
                status_dict[dep_id] == State.failed
                for dep_id in request.deps
            )

            if any_dep_failed:
                status_dict[task_id] = State.failed
                logger.error(
                    'Task %s failed since a dependency failed.',
                    task_id,
                    exc_info=True
                )
                continue

            has_non_satisfied_dep = False
            for dep_id in request.deps:
                if status_dict[dep_id] != State.completed:
                    logger.debug(
                        'Task %s set to wait for %s.',
                        task_id,
                        dep_id,
                    )

                    if dep_id not in deps_dict:
                        deps_dict[dep_id] = []

                    # Do this weird thing to sync the dictionary. Using
                    # deps_dict[dep_id].append(...) doesn't work.
                    tmp = deps_dict[dep_id]
                    tmp.append((task_id, request))
                    deps_dict[dep_id] = tmp

                    has_non_satisfied_dep = True

            if has_non_satisfied_dep:
                continue

        logger.debug(
            'Task %s started target %r.',
            task_id, request.target
        )

        status_dict[task_id] = State.started
        try:
            process = subprocess.Popen(
                ['bash'],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                cwd=request.target.working_dir,
            )

            stdout, stderr = process.communicate(request.target.spec)

            with open(_get_stdout_path(request.target), 'w') as stdout_fileobj:
                stdout_fileobj.write(stdout)
            with open(_get_stderr_path(request.target), 'w') as stderr_fileobj:
                stderr_fileobj.write(stderr)

            if process.returncode != 0:
                raise Exception(stderr)
        except:
            status_dict[task_id] = State.failed

            logger.error(
                'Task %s failed.',
                task_id,
                exc_info=True
            )
        else:
            status_dict[task_id] = State.completed
            logger.debug(
                'Task %s completed target %r.',
                task_id, request.target
            )
        finally:
            with deps_lock:
                if task_id in deps_dict:
                    logger.debug(
                        'Task %s has waiting dependents. Requeueing.', task_id)
                    for dep_task_id, dep_request in deps_dict[task_id]:
                        task_queue.put((dep_task_id, dep_request))


def handle_request(request, task_queue, status_dict):
    try:
        logger.debug('Received request %r.', request)
        return request.handle(task_queue, status_dict)
    except:
        logger.error('Invalid request %r.', request, exc_info=True)


def handle_client(conn, task_queue, status_dict):
    logger.debug('Accepted client connection.')
    try:
        while True:
            request = conn.recv()
            response = handle_request(request, task_queue, status_dict)
            if response is not None:
                conn.send(response)
    except EOFError:
        logger.debug('Client connection closed.')


def wait_for_clients(address, task_queue, status_dict):
    serv = Listener(address)
    while True:
        try:
            client = serv.accept()
            handle_client(client, task_queue, status_dict)
        except Exception:
            traceback.print_exc()


def start(hostname='', port=25000, num_workers=None):
    try:
        with Manager() as manager:
            status_dict = manager.dict()
            deps_dict = manager.dict()
            task_queue = manager.Queue()
            deps_lock = manager.Lock()

            workers = Pool(
                processes=num_workers,
                initializer=worker,
                initargs=(task_queue, status_dict, deps_dict, deps_lock),
            )

            wait_for_clients((hostname, port), task_queue, status_dict)
    except KeyboardInterrupt:
        logger.debug('Shutting down...')
        workers.close()
        sys.exit(0)

if __name__ == '__main__':
    start(num_workers=8)
