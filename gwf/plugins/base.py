from ..ext import Extension


class Plugin(Extension):
    """Abstract base class for plugins."""

    def configure(self, workflow, backend, config):
        """Called to configure the plugin.

        This method sets the workflow, backend and config attributes on the
        object, but may be overridden to perform other kinds of initialization.
        """
        self.workflow = workflow
        self.backend = backend
        self.config = config
