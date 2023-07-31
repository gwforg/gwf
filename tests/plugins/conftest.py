import pytest

from gwf.backends.local import start_cluster_in_background
from gwf.conf import FileConfig


@pytest.fixture
def local_backend(tmp_path):
    config = FileConfig.load(tmp_path.joinpath(".gwfconf.json"))
    config["backend"] = "local"
    config.dump()

    with start_cluster_in_background(tmp_path, 1, "localhost", 12345):
        yield
