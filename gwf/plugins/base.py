from ..ext import Extension


class Plugin(Extension):
    """Abstract base class for plugins."""

    def configure(self, get_graph, get_active_backend, config):
        """Called to configure the plugin.

        This method sets the workflow, backend and config attributes on the
        object, but may be overridden to perform other kinds of initialization.
        """
        self.get_graph = get_graph
        self.get_active_backend = get_active_backend
        self.config = config
