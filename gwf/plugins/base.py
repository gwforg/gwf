from ..ext import Extension


class Plugin(Extension):
    """Abstract base class for plugins."""

    def configure(self, get_prepared_workflow, get_active_backend, config):
        """Called to configure the plugin.

        This method sets the workflow, backend and config attributes on the
        object, but may be overridden to perform other kinds of initialization.
        """
        self.get_prepared_workflow = get_prepared_workflow
        self.get_active_backend = get_active_backend
        self.config = config
