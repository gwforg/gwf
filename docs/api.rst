API
===

The implementation of *gwf* consists of a few main abstractions. Units of work
are defined by creating :class:`~gwf.Target` instances which also define
the files used and produced by the target. A :class:`~gwf.Workflow` ties
together and allows for easy creation of targets.

When all targets have been defined on a workflow, the workflow must be turned
into a :class:`~gwf.PreparedWorkflow`. Preparing the workflow consists of
computing the entire dependency graph of the workflow, checking the workflow
for inconsistensies and circular dependencies, and figuring out which targets
should be run (see :meth:`~gwf.PreparedWorkflow.should_run`).

A :class:`~gwf.PreparedWorkflow` can be given to a :class:`~gwf.Backend`
implementation which will then execute the workflow.

Core
----

.. automodule:: gwf
   :members:

Backends
--------

.. autoclass:: gwf.backends.base.Backend
   :members:
