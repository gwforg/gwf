from ..ext import Extension


class Plugin(Extension):
    """Abstract base class for plugins."""

    def configure(self, workflow, backend, config):
        self.workflow = workflow
        self.backend = backend
        self.config = config
