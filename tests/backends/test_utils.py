import pytest

from gwf.backends.exceptions import BackendError
from gwf.backends.utils import call


def test_call_with_missing_executable_names_it_without_assuming_slurm(monkeypatch):
    monkeypatch.setenv("PATH", "")

    with pytest.raises(BackendError) as excinfo:
        call("qsub", "-terse")

    message = str(excinfo.value)
    assert '"qsub"' in message
    assert "Slurm" not in message
