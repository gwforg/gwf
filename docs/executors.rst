.. _executors:

Executors
=========

.. note::

   Executors are currently only available for the Slurm backend. Setting an 
   executor for other backends will not have any effect. Other backends will
   be supported in the near future.

Executors are used to enable runtime behavior for targets. This means that you
can now use executors to run your target inside Conda, Pixi, Apptainer, and
Singularity environments:

.. code-block:: python

   from gwf import Workflow
   from gwf.executors import Conda

   gwf = Workflow()

   gwf.target("Test", inputs=[], outputs=[], executor=Conda("myenv")) <<< """
   echo this will run inside the `myenv` Conda environment.
   """

A default executor can also be specified for the workflow:

.. code-block:: python

   from gwf import Workflow
   from gwf.executors import Conda

   gwf = Workflow(executor=Conda("myenv"))

   gwf.target("Test", inputs=[], outputs=[]) <<< """
   echo this will run inside the `myenv` Conda environment.
   """

Available executors
-------------------

.. automodule:: gwf.executors
   :members:
   :exclude-members: Executor
   :member-order: bysource
