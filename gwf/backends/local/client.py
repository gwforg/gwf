from multiprocessing.connection import Client as Client_

from .server import StatusRequest, SubmitRequest


class Client:

    def __init__(self, *args, **kwargs):
        self.client = Client_(*args, **kwargs)

    def submit(self, target, deps=None):
        if deps is None:
            deps = []
        request = SubmitRequest(target=target, deps=deps)
        self.client.send(request)
        return self.client.recv()

    def status(self):
        request = StatusRequest()
        self.client.send(request)
        return self.client.recv()

    def close(self):
        self.client.close()


__all__ = ('Client',)
