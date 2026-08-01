"""gwf - a pragmatic workflow tool"""

from .core import AnonymousTarget, Target
from .exec import IsolationConfig
from .workflow import TargetList, Workflow

__version__ = "2.1.1"

__all__ = (
    "Target",
    "AnonymousTarget",
    "Workflow",
    "TargetList",
    "IsolationConfig",
)
