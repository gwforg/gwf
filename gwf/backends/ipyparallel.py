import json
import logging
import os
import os.path
import subprocess
import warnings

from . import Backend
from ..exceptions import BackendError

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from ipyparallel import Client


warnings.simplefilter("ignore")

logger = logging.getLogger(__name__)


def _get_log_dir_path(target):
    return os.path.join(target.working_dir, '.gwf', 'logs')


def _get_stdout_path(target):
    return os.path.join(_get_log_dir_path(target), '{}.stdout'.format(target.name))


def _get_stderr_path(target):
    return os.path.join(_get_log_dir_path(target), '{}.stderr'.format(target.name))


def _make_log_dir(target):
    log_dir = _get_log_dir_path(target)
    try:
        os.makedirs(log_dir)
    except OSError:
        pass
    return log_dir


def _run_target(target):
    process = subprocess.Popen(
        ['bash'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        cwd=target.working_dir,
    )

    stdout, stderr = process.communicate(target.spec)

    with open(_get_stdout_path(target), 'w') as stdout_fileobj:
        stdout_fileobj.write(stdout)

    with open(_get_stderr_path(target), 'w') as stderr_fileobj:
        stderr_fileobj.write(stderr)

    if process.returncode != 0:
        raise BackendError(stderr)


class IPyParallelBackend(Backend):
    """A backend that runs targets on an iPyParallel cluster."""

    name = 'ipyparallel'

    supported_options = []
    option_defaults = {}

    def configure(self, *args, **kwargs):
        super().configure(*args, **kwargs)
        self.client = Client()
        self.view = self.client.load_balanced_view()
        self.results = {}

        try:
            with open('.gwf/ipyparallel-db.json', 'r') as fileobj:
                self.results = json.load(fileobj)
        except OSError:
            logger.debug(
                'Could not open job database. Assuming it is empty.'
            )
            self.results = {}

        self.results = {
            k: v
            for k, v in self.results.items()
            if k in self.workflow.targets and v in self.client.outstanding
        }

    def submit(self, target):
        deps = [
            self.results[dep]
            for dep in self.workflow.dependencies[target]
            if dep in self.results
        ]

        with self.view.temp_flags(after=deps, block=False):
            result = self.view.apply(_run_target, target)

            # Our messages never have any children, so we always only have one
            # message id.
            msg_id = result.msg_ids[0]
            logger.debug(
                'Target "{}" submitted, received id {}.'.format(
                    target.name, msg_id,
                )
            )
            self.results[target.name] = msg_id

    def cancel(self, target):
        msg_id = self.results.get(target.name)
        if msg_id is None:
            return
        self.client.abort(jobs=[msg_id])

    def submitted(self, target):
        return target.name in self.results

    def running(self, target):
        msg_id = self.results.get(target.name)
        if msg_id is None:
            return False
        return (msg_id in self.client.outstanding and
                msg_id not in self.client.results)

    def logs(self, target, stderr=False):
        if stderr:
            stderr_path = _get_stderr_path(target)
            return open(stderr_path)
        stdout_path = _get_stdout_path(target)
        return open(stdout_path)

    def close(self):
        self.client.close()
        with open('.gwf/ipyparallel-db.json', 'w') as fileobj:
            json.dump(self.results, fileobj)
