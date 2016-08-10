BACKENDS = {}


def register_backend(backend):
    BACKENDS[backend.name] = backend


class BackendType(type):

    def __new__(meta, name, bases, class_dict):
        cls = type.__new__(meta, name, bases, class_dict)
        register_backend(cls)
        return cls


class Backend(metaclass=BackendType):

    def submitted(self, targets):
        raise NotImplementedError()

    def running(self, targets):
        raise NotImplementedError()

    def submit(self, targets):
        raise NotImplementedError()

    def cancel(self, targets):
        raise NotImplementedError()
