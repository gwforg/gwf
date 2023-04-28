import threading
import pytest
from gwf.backends.local import Cluster


@pytest.fixture
def local_backend():
    cluster = Cluster(num_workers=1)
    thread = threading.Thread(target=cluster.start)
    thread.start()
    yield
    cluster.shutdown()
    thread.join()