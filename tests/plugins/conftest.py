import tempfile
import time
from pathlib import Path

import pytest

from gwf.conf import FileConfig
from tests.local_backend import background_scheduler


@pytest.fixture
def local_backend(tmp_path, monkeypatch):
    with tempfile.TemporaryDirectory(prefix="gwf-test-", dir="/tmp") as runtime_dir:
        monkeypatch.setenv("XDG_RUNTIME_DIR", runtime_dir)

        config = FileConfig.load(tmp_path.joinpath(".gwfconf.json"))
        config["backend"] = "local"
        config.dump()

        with background_scheduler(tmp_path, 1):
            yield


@pytest.fixture
def wait_for_path():
    def wait(path, timeout=10):
        path = Path(path)
        deadline = time.monotonic() + timeout
        while not path.exists():
            if time.monotonic() >= deadline:
                pytest.fail(f"Timed out waiting for {path}")
            time.sleep(0.1)

    return wait
