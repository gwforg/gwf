import tempfile

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
