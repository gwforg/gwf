import shutil
import subprocess

from .exceptions import BackendError


def _find_exe(name):
    exe = shutil.which(name)
    if exe is None:
        raise BackendError(
            f'Could not find executable "{name}". Make sure that the workload '
            f"manager used by this backend is installed on this host and that "
            f'"{name}" is on your PATH.'
        )
    return exe


def has_exe(name):
    return shutil.which(name) is not None


def call(executable_name, *args, input=None, environ=None):
    executable_path = _find_exe(executable_name)
    proc = subprocess.Popen(
        [executable_path] + list(args),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,
        universal_newlines=True,
        env=environ,
    )
    stdout, stderr = proc.communicate(input)

    # Some commands, like scancel, do not return a non-zero exit code if they
    # fail. The only way to check if they failed is by checking whether an
    # error message occurred in standard error, so we check both the return
    # code and stderr.
    if proc.returncode != 0 or "error:" in stderr:
        raise BackendError(stderr)
    return stdout
