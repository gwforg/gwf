.. _tutorial:

Tutorial
========

In this tutorial we will explore various concepts in *gwf*. We will define
workflows and see how *gwf* can help us keep track of the progress of workflow
execution, the output of targets and dependencies between targets. Have fun!

We'll assume that you have the Anaconda_ distribution installed and that you are
familiar with how to install and manage packages and environments through the
*conda* package manager.

First, let's install *gwf* in its own conda environment. Create a new environment
for your project, we'll call it *myproject*.

.. code-block:: console

    $ conda create -n myproject python=3.4
    $ source activate myproject
    $ conda install -c dansondergaard gwf

You should now be able to run the following command.

.. code-block:: console

    $ gwf -h

This should show you the commands and options available through *gwf*. If you just
run:

.. code-block:: console

    $ gwf

you'll get an error that looks something like this:

.. code-block:: console

    ERROR |  The file "/Users/das/Desktop/test-gwf/workflow.py" does not exist.

We get this error since we didn't define a workflow file yet.

A Minimal Workflow
------------------

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

There you go! We have now declared a workflow with one target and that target
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
Then run the command ``gwf -b local workers``. This will start a pool of workers that
*gwf* can now submit targets to.

Switch back to the other terminal and then run:

.. code-block:: console

    $ gwf -b local run

*gwf* schedules and then submits ``MyTarget`` to the pool of workers you started in
the other terminal window (the ``-b local`` flag tells *gwf* to use the
:class:`~gwf.backends.local.LocalBackend`). This command doesn't output anything
since *gwf* tries to only output something when explicitly told so, or if something
is wrong.

Within a few seconds you should see ``greeting.txt`` in the project directory. Try
to open it in your favorite text editor!

To actually see what happens when you run ``gwf -b local run``, try to delete
``greeting.txt`` and then run:

.. code-block:: console

    $ gwf -b local -v info run

The ``-v info`` flag tells *gwf* to output a bit more information when it runs.
If you want even more information you may use ``-v debug``. The command show now
output this:

.. code-block:: console

    INFO  |  Scheduling target MyTarget.
    INFO  |  Submitting target MyTarget.

This says that *gwf* considered the target for execution and then decided to submit
it to the backend (in this case because the output file, ``greeting.txt``, does not
already exist). After a few seconds, you should see that ``greeting.txt`` has been
created again.

Now try the same command again:

.. code-block:: console

    $ gwf -b local -v info run
    INFO  |  Scheduling target MyTarget.

This time, *gwf* considers the target for submission, but decides not to submit it
since all of the output files (only one in this case) exist.

Setting a Default Backend
-------------------------

By now you probably got really tired of typing ``-b local`` for every single command.
To save some keystrokes, let's set the local backend as the default for this project.

.. code-block:: console

    $ gwf config backend local

This creates a configuration file in the project folder and sets the backend to
``local`` by default. To test it out, let's try to run the same command as before,
but without the ``-b local`` flag.

.. code-block:: console

    $ gwf -v info run
    INFO  |  Scheduling target MyTarget.

*gwf* now uses the local backend by default, so everything works as before. If you
are crazy about seeing what *gwf* does, you can also get rid of the ``-v info``
flag by setting the default verbosity level.

.. code-block:: console

    $ gwf config verbosity info
    $ gwf run
    INFO  |  Scheduling target MyTarget.

As we'd expect, *gwf* outputs the same as before, but this time we didn't have to
set the ``-v info`` flag!

Defining Targets With Dependencies
----------------------------------



Observing Targets Execution
---------------------------

Status, logs.


Reusable Targets With Templates
-------------------------------


Cleaning Up
-----------


Running with Another Backend
----------------------------


.. _Anaconda: https://www.continuum.io/downloads
