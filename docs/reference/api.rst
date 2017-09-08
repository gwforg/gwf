API
===

The implementation of *gwf* consists of a few main abstractions. Units of work
are defined by creating :class:`~gwf.Target` instances which also define
the files used and produced by the target. A :class:`~gwf.Workflow` ties
together and allows for easy creation of targets.

When all targets have been defined on a workflow, the workflow is turned
into a :class:`~gwf.Graph` which will compute the entire dependency graph of the
workflow, checking the workflow for inconsistencies and circular dependencies.

A target in a :class:`~gwf.Graph` can be scheduled on a
:class:`~gwf.backends.Backend` using the :class:`~gwf.Scheduler`.

.. automodule:: gwf
   :members: Target, Workflow, Graph, Scheduler

.. autoclass:: gwf.backends.Backend
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: gwf.backends.Status
   :members:
   :undoc-members:
