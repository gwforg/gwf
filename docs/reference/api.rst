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

Workflow
--------

.. automodule:: gwf
    :members: Workflow, Target, AnonymousTarget, TargetList

Core
----

.. automodule:: gwf.core
    :members: Graph, Scheduler

Backends
--------

.. automodule:: gwf.backends
    :members: Backend, Status
    :inherited-members:
    :undoc-members:


Log Managers
~~~~~~~~~~~~

.. automodule:: gwf.backends.logmanager
    :members: FileLogManager, MemoryLogManager

Filtering
---------

.. automodule:: gwf.filtering
    :members: filter_generic, filter_names
    :undoc-members:

Helpers for filtering:

.. automodule:: gwf.filtering
    :members: ApplyMixin
    :undoc-members:
