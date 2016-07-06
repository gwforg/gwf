import gwf

from gwf.api import target, template, memorize, include_workflow, shell, glob


__all__ = ['target', 'template', 'memorize', 'include_workflow', 'shell', 'glob']

# This global variable will hold all the targets after a workflow script has completed.
# gwf will use this list for its further processing.
ALL_TARGETS = {}

# This will be set in the gwf script and refer to the grid backend used.
BACKEND = None

