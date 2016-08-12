BACKENDS = {}


def register_backend(name, backend_cls):
    BACKENDS[name] = backend_cls


class BackendType(type):

    def __new__(meta, name, bases, class_dict):
        cls = type.__new__(meta, name, bases, class_dict)
        register_backend(name, cls)
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
