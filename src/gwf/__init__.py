"""gwf - a pragmatic workflow tool"""

from .core import AnonymousTarget, Module, Target
from .temp import temp
from .workflow import TargetList, Workflow

__version__ = "2.1.1"

__all__ = ("Target", "AnonymousTarget", "Module", "Workflow", "TargetList", "temp")
