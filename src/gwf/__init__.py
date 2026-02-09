"""gwf - a pragmatic workflow tool"""

from .core import AnonymousTarget, Target
from .path import protect, temp
from .workflow import TargetList, Workflow

__version__ = "2.1.1"

__all__ = (
    "Target",
    "AnonymousTarget",
    "Workflow",
    "TargetList",
    "protect",
    "temp",
)
