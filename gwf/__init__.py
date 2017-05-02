from .core import Graph, Target, Workflow

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


__all__ = ('Target', 'Workflow', 'Graph',)
