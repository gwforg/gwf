API
===

The implementation of *gwf* consists of a few main abstractions. Units of work
are defined by creating :class:`~gwf.Target` instances which also define
the files used and produced by the target. A :class:`~gwf.Workflow` ties
together and allows for easy creation of targets.

When all targets have been defined on a workflow, the workflow must be turned
into a :class:`~gwf.Graph`. Preparing the workflow consists of
computing the entire dependency graph of the workflow, checking the workflow
for inconsistencies and circular dependencies.

A targets in a :class:`~gwf.Graph` can be scheduled on a
:class:`~gwf.backends.Backend` using the :func:`gwf.core.schedule_many`
function.


Core
----

.. automodule:: gwf
   :members: Target, Workflow, Graph, schedule_many


Backends
--------

.. autoclass:: gwf.backends.base.Backend
   :members:
   :inherited-members:
   :undoc-members:


Exceptions
----------

.. automodule:: gwf.exceptions
   :members:
   :member-order: bysource
   :undoc-members:
