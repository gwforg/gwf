import subprocess
from distutils.spawn import find_executable

from .exceptions import BackendError


def _find_exe(name):
    exe = find_executable(name)
    if exe is None:
        raise BackendError(
            'Could not find executable "{}". This backend requires Slurm to be installed on this host.'.format(
                name
            )
        )
    return exe


def call(executable_name, *args, input=None):
    executable_path = _find_exe(executable_name)
    proc = subprocess.Popen(
        [executable_path] + list(args),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = proc.communicate(input)

    # Some commands, like scancel, do not return a non-zero exit code if they
    # fail. The only way to check if they failed is by checking whether an
    # error message occurred in standard error, so we check both the return
    # code and stderr.
    if proc.returncode != 0 or "error:" in stderr:
        raise BackendError(stderr)
    return stdout
