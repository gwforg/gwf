.. _tutorial:

Tutorial
========

In this tutorial we will explore various concepts in *gwf*. We will define
workflows and see how *gwf* can help us keep track of the progress of workflow
execution, the output of targets and dependencies between targets. Have fun!

A Minimal Example
-----------------

To get started we must define a *workflow file* containing a workflow to which
we can add targets. Unless *gwf* is told otherwise it assumes that the workflow
file is called ``workflow.py`` and that the workflow is called `gwf` (see
`Defining Multiple Workflows in a File`_ for an example where a workflow is not
necessarily called `gwf`)::

    from gwf import Workflow

    gwf = Workflow()

    gwf.target('MyTarget') << """
    echo hello world
    """

In the example above we define a workflow and then add a target called
``MyTarget``. A target is a single unit of computation that uses zero or more
files (inputs) and produces zero or more files (outputs).

The target defined above does not use any files and doesn't produce any files
either. However, it does run a single command (``echo hello world``), but the
output of the command is thrown away. Let's fix that! Change the target
definition to this::

    gwf.target('MyTarget', outputs=['greeting.txt']) << """
    echo hello world
    """

This tells *gwf* that the target will create a file called ``greeting.txt`` when
it is run. However, the target does not actually create the file yet. Let's
fix that too::

    gwf.target('MyTarget', outputs=['greeting.txt']) << """
    echo hello world > greeting.txt
    """

There ya' go! We have now declared a workflow with one target and that target
creates the file ``greeting.txt`` with the line ``hello world`` in it. Now let's
try to run our workflow...

Running Your First Workflow
---------------------------

First, let's make a directory for our project. We'll call the directory
``myproject``. Now create an empty file called ``workflow.py`` in the project
directory and paste the workflow specification into it::

    from gwf import Workflow

    gwf = Workflow()

    gwf.target('MyTarget', outputs=['greeting.txt']) << """
    echo hello world > greeting.txt
    """

We're now ready to run our workflow. However, *gwf* does not actually execute
the targets in a workflow, it only schedules the target using a backend. This
may sound cumbersome, but it enables *gwf* to run workflows in very different
environments: anything from your laptop to a cluster with thousands of cores
available.

For this tutorial we just want to run our workflows locally. To do this we can
use the built-in ``local`` backend. Essentially this backend allows you to run
workflows utilizing all cores of your computer and thus it can be very useful
for small workflows that don't require a lot of resources.

First, open another terminal window and navigate to the ``myproject`` directory.
Then run the command ``gwf workers``. This will start a pool of workers that
*gwf* can now submit targets to.

Switch back to the other terminal and then run the command ``gwf run``. If
everything is fine, you should see output like this::

    ### EXAMPLE OUTPUT ###

We can see that *gwf* scheduled and then submitted ``MyTarget``.
Within a few seconds you should see ``greeting.txt`` in the project directory.

What happens if we type ``gwf run`` again? Let's try!::

    ### EXAMPLE OUTPUT FROM SECOND RUN ###

This time the target is not submitted to the backend since ``greeting.txt``
already exists and is up to date. Try to delete ``greeting.txt``, then run
the workflow again::

    ## EXAMPLE OUTPUT FROM SECOND RUN ###

The target was submitted again, what a relief!

Observing Running Targets
-------------------------

Status, logs.


Defining Targets With Dependencies
----------------------------------


Reusable Targets With Templates
-------------------------------


Defining Multiple Workflows in a File
-------------------------------------


Including Other Workflows
-------------------------


Running with Another Backend
----------------------------
