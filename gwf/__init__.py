from .core import Graph, Target, Workflow, schedule_many

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


__all__ = ('Target', 'Workflow', 'Graph', 'schedule_many',)
