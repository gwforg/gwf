from unittest import mock
from unittest.mock import call

import pytest

from gwf import Target
from gwf.backends import Status
from gwf.backends.exceptions import BackendError, DependencyError, TargetError
from gwf.backends.slurm import SlurmBackend
from gwf.conf import config


@pytest.fixture(autouse=True)
def setup(monkeypatch, tmpdir):
    exes = {"sbatch": "/bin/sbatch", "squeue": "/bin/squeue", "scancel": "/bin/scancel"}
    monkeypatch.setattr("gwf.backends.utils.find_executable", lambda e: exes[e])


@pytest.fixture
def popen(request):
    patched = mock.patch("gwf.backends.utils.subprocess.Popen")
    mocked = patched.__enter__()

    @request.addfinalizer
    def unpatch():
        patched.__exit__()

    return mocked


def test_initialization(tmpdir, popen, no_sleep):
    t1 = Target.empty("Target1")
    t2 = Target.empty("Target2")
    t3 = Target.empty("Target3")
    t4 = Target.empty("Target4")

    gwf_dir = tmpdir.ensure_dir(".gwf")
    tracked_file = gwf_dir.join("slurm-backend-tracked.json")
    tracked_file.write(
        '{"Target1": "1", "Target2": "2", "Target3": "3", "Target4": "4"}'
    )

    popen.return_value.returncode = 0
    popen.return_value.communicate.return_value = ("3;SE\n2;BF\n1;R\n", "")

    with tmpdir.as_cwd():
        backend = SlurmBackend()
        assert backend.status(t1) == Status.RUNNING
        assert backend.status(t2) == Status.UNKNOWN
        assert backend.status(t3) == Status.SUBMITTED
        assert backend.status(t4) == Status.UNKNOWN

        popen.return_value.returncode = 1
        popen.return_value.communicate.return_value = ("", "This is stderr")
        with pytest.raises(BackendError):
            with SlurmBackend():
                pass


def test_submit_1(popen, no_sleep):
    t1 = Target.empty("Target1")
    t2 = Target.empty("Target2")
    t3 = Target.empty("Target3")
    t4 = Target.empty("Target4")

    popen.return_value.returncode = 0
    popen.return_value.communicate.return_value = ("", "")
    backend = SlurmBackend()

    popen.return_value.returncode = 0
    popen.return_value.communicate.side_effect = [
        ("1\n", ""),
        ("2\n", ""),
        ("3\n", ""),
        ("", "This is stderr\n"),
    ]
    backend.submit(t1, [])
    backend.submit(t2, [])
    backend.submit(t3, [t1, t2])

    assert (
        call(
            ["/bin/sbatch", "--parsable"],
            stderr=-1,
            stdin=-1,
            stdout=-1,
            universal_newlines=True,
        )
        in popen.call_args_list
    )

    assert (
        call(
            ["/bin/sbatch", "--parsable"],
            stderr=-1,
            stdin=-1,
            stdout=-1,
            universal_newlines=True,
        )
        in popen.call_args_list
    )

    assert (
        call(
            ["/bin/sbatch", "--parsable", "--dependency=afterok:1:2"],
            stderr=-1,
            stdin=-1,
            stdout=-1,
            universal_newlines=True,
        )
        in popen.call_args_list
    )

    with pytest.raises(DependencyError):
        backend.submit(t3, [t4])

    popen.return_value.returncode = 1
    popen.return_value.communicate.side_effect = [
        ("", "This is stderr\n"),
        ("", "This is stderr\n"),
        ("", "This is stderr\n"),
    ]
    with pytest.raises(BackendError):
        backend.submit(t4, [])


def test_submit_2(popen, monkeypatch, no_sleep):
    t = Target(
        name="TestTarget",
        inputs=[],
        outputs=[],
        working_dir="/some/dir",
        options={
            "cores": 16,
            "memory": "16g",
            "walltime": "12:00:00",
            "queue": "normal",
            "account": "someaccount",
            "constraint": "graphics*4",
            "mail_type": "BEGIN,END,FAIL",
            "mail_user": "test@domain.com",
            "qos": "somename",
        },
        spec="echo hello world",
    )

    popen.return_value.returncode = 0
    popen.return_value.communicate.return_value = ("", "")
    backend = SlurmBackend()

    backend.submit(t, [])
    (script,), _ = popen.return_value.communicate.call_args

    assert "#!/bin/bash" in script
    assert "#SBATCH --job-name=TestTarget" in script
    assert "#SBATCH -c 16" in script
    assert "#SBATCH --mem=16g" in script
    assert "#SBATCH -t 12:00:00" in script
    assert "#SBATCH -p normal" in script
    assert "#SBATCH -A someaccount" in script
    assert "#SBATCH -C graphics*4" in script
    assert "#SBATCH --mail-type=BEGIN,END,FAIL" in script
    assert "#SBATCH --mail-user=test@domain.com" in script
    assert "#SBATCH --qos=somename" in script
    assert "#SBATCH --output=.gwf/logs/TestTarget.stdout" in script
    assert "#SBATCH --error=.gwf/logs/TestTarget.stderr" in script
    assert "cd /some/dir" in script
    assert "export GWF_JOBID=$SLURM_JOBID" in script
    assert "echo hello world" in script

    monkeypatch.setitem(config, "backend.slurm.log_mode", "merged")

    backend.submit(t, [])
    (script,), _ = popen.return_value.communicate.call_args

    assert "#SBATCH --output=.gwf/logs/TestTarget.stdout" in script
    assert "#SBATCH --error=.gwf/logs/TestTarget.stderr" not in script

    monkeypatch.setitem(config, "backend.slurm.log_mode", "none")

    backend.submit(t, [])
    (script,), _ = popen.return_value.communicate.call_args

    assert "#SBATCH --output=/dev/null" in script
    assert "#SBATCH --error=" not in script


def test_cancel(popen, no_sleep):
    t = Target.empty("Target")

    popen.return_value.returncode = 0
    popen.return_value.communicate.return_value = ("", "")
    backend = SlurmBackend()
    with pytest.raises(TargetError):
        backend.cancel(t)

    popen.return_value.communicate.return_value = ("1\n", "")
    backend.submit(t, [])
    assert backend.status(t) == Status.SUBMITTED

    backend.cancel(t)
    assert backend.status(t) == Status.UNKNOWN

    assert (
        call(
            ["/bin/scancel", "--verbose", "1"],
            stderr=-1,
            stdin=-1,
            stdout=-1,
            universal_newlines=True,
        )
        in popen.call_args_list
    )

    popen.return_value.returncode = 1
    popen.return_value.communicate.return_value = ("", "This is stderr")
    with pytest.raises(BackendError):
        backend.cancel(t)


def test_close(tmpdir, popen, no_sleep):
    tmpdir.ensure_dir(".gwf")
    with tmpdir.as_cwd():
        popen.return_value.returncode = 0
        popen.return_value.communicate.return_value = ("", "")

        backend = SlurmBackend()
        backend.close()
    assert tmpdir.join(".gwf", "slurm-backend-tracked.json").check()
