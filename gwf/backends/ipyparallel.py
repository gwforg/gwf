import logging
import os
import os.path
import subprocess

from ipyparallel import Client

from . import Backend

logger = logging.getLogger(__name__)


def _make_log_dir(target):
    log_dir = os.path.join(target.working_dir, '.gwf', 'logs')
    try:
        os.makedirs(log_dir)
    except OSError:
        pass
    return log_dir


def _run_target(target):
    log_dir = _make_log_dir(target)

    process = subprocess.Popen(
        ['bash'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        cwd=target.working_dir,
    )

    stdout, stderr = process.communicate(target.spec)

    stdout_filename = '{}.stdout'.format(target.name)
    stdout_path = os.path.join(log_dir, stdout_filename)
    with open(stdout_path, 'w') as stdout_fileobj:
        stdout_fileobj.write(stdout)

    stderr_filename = '{}.stderr'.format(target.name)
    stderr_path = os.path.join(log_dir, stderr_filename)
    with open(stderr_path, 'w') as stderr_fileobj:
        stderr_fileobj.write(stderr)

    if process.returncode != 0:
        raise Exception(stderr)


class IPyParallelBackend(Backend):

    name = 'ipyparallel'
    supported_options = []

    def configure(self, *args, **kwargs):
        super().configure(*args, **kwargs)
        self.client = Client()
        self.view = self.client.load_balanced_view()
        self.results = {}

    def submit(self, target):
        logger.debug(target)
        logger.debug('self.results = %r', self.results)
        deps = []
        for dep in self.workflow.dependencies[target]:
            if dep in self.results:
                deps.append(self.results[dep])

        with self.view.temp_flags(after=deps, block=False):
            self.results[target] = self.view.apply(_run_target, target)

    def cancel(self, target):
        pass

    def submitted(self, target):
        return False

    def running(self, target):
        return False

    def logs(self, target, stderr=False):
        return ''

    def close(self):
        print(self.view.wait(self.results.values()))
