from contextlib import contextmanager
import multiprocessing

from gwf.backends.exceptions import BackendError
from gwf.backends.local import Client, get_socket_path, start_cluster


def _stop_process(process):
    process.join(timeout=5)
    if process.is_alive():
        process.kill()
        process.join()


@contextmanager
def background_scheduler(working_dir, max_cores, **kwargs):
    socket_path = get_socket_path(working_dir)
    ready_event = multiprocessing.Event()
    kwargs["ready_event"] = ready_event
    process = multiprocessing.Process(
        target=start_cluster,
        args=(working_dir, max_cores),
        kwargs=kwargs,
    )
    process.start()

    if not ready_event.wait(timeout=5):
        _stop_process(process)
        raise BackendError("Local scheduler did not start.")

    try:
        yield socket_path
    finally:
        try:
            if process.is_alive():
                with Client.connect(socket_path) as client:
                    client.shutdown()
        finally:
            _stop_process(process)
