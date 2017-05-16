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

    $ gwf

This should show you the commands and options available through *gwf*.

.. caution::
    You may see an error similar to this when you try running *gwf*::

        UnicodeEncodeError: 'charmap' codec can't encode character '\u2020' in position 477: character maps to <undefined>

    This error occurs because your isn't configured to use UTF-8 as the default encoding.
    To fix the error insert the following lines in your ``.bashrc`` file::

        export LANG=en_US.utf8
        export LC_ALL=en_US.utf8

    If you're not in the US you may want to set it to something else. For example,
    if you're in Denmark you may want to use the following configuration::

        export LANG=da_DK.utf8
        export LC_ALL=da_DK.utf8


A Minimal Workflow
------------------

To get started we must define a *workflow file* containing a workflow to which
we can add targets. Unless *gwf* is told otherwise it assumes that the workflow
file is called ``workflow.py`` and that the workflow is called `gwf`::

    from gwf import Workflow

    gwf = Workflow()

    gwf.target('MyTarget', inputs=[], outputs=[]) << """
    echo hello world
    """

In the example above we define a workflow and then add a target called
``MyTarget``. A target is a single unit of computation that uses zero or more
files (inputs) and produces zero or more files (outputs).

The target defined above does not use any files and doesn't produce any files
either. However, it does run a single command (``echo hello world``), but the
output of the command is thrown away. Let's fix that! Change the target
definition to this::

    gwf.target('MyTarget', inputs=[], outputs=['greeting.txt']) << """
    echo hello world
    """

This tells *gwf* that the target will create a file called ``greeting.txt`` when
it is run. However, the target does not actually create the file yet. Let's
fix that too::

    gwf.target('MyTarget', inputs=[], outputs=['greeting.txt']) << """
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

    gwf.target('MyTarget', inputs=[], outputs=['greeting.txt']) << """
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
Then run the command:

.. code-block:: console

    $ gwf -b local workers

This will start a pool of workers that *gwf* can now submit targets to.
Switch back to the other terminal and then run:

.. code-block:: console

    $ gwf run

*gwf* schedules and then submits ``MyTarget`` to the pool of workers you started in
the other terminal window (the ``-b local`` flag tells *gwf* to use the
:class:`~gwf.backends.local.LocalBackend`). This command doesn't output anything
since *gwf* tries to only output something when explicitly told so, or if something
is wrong.

Within a few seconds you should see ``greeting.txt`` in the project directory. Try
to open it in your favorite text editor!

To actually see what happens when you run ``gwf run``, try to delete
``greeting.txt`` and then run:

.. code-block:: console

    $ gwf -v info run

The ``-v info`` flag tells *gwf* to output a bit more information when it runs.
If you want even more information you may use ``-v debug``. The command should now
output this:

.. code-block:: console

    Scheduling target MyTarget.
    Submitting target MyTarget.

This says that *gwf* considered the target for execution and then decided to submit
it to the backend (in this case because the output file, ``greeting.txt``, does not
already exist). After a few seconds, you should see that ``greeting.txt`` has been
created again.

Now try the same command again:

.. code-block:: console

    $ gwf -b local -v info run
    Scheduling target MyTarget.

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
    Scheduling target MyTarget.

*gwf* now uses the local backend by default, so everything works as before. If you
are crazy about seeing what *gwf* does, you can also get rid of the ``-v info``
flag by setting the default verbosity level.

.. code-block:: console

    $ gwf config verbosity info
    $ gwf run
    Scheduling target MyTarget.

As we'd expect, *gwf* outputs the same as before, but this time we didn't have to
set the ``-v info`` flag!

From now on we'll assume that you've started a pool of workers for the local backend
and configured this to be the default backend.

Defining Targets with Dependencies
----------------------------------

Targets in *gwf* represent isolated units of work. However, we can declare
dependencies between targets to construct complex workflows. A target B that
depends on a target A will only run when A has been run successfully (that
is, if all of the output files of A exist).

In *gwf*, dependencies are declared through file dependencies. This is best
understood through an example::

    from gwf import Workflow

    gwf = Workflow()

    gwf.target('TargetA', inputs=[], outputs=['x.txt']) << """
    echo "this is x" > x.txt
    """

    gwf.target('TargetB', inputs=[], outputs=['y.txt']) << """
    echo "this is y" > y.txt
    """

    gwf.target('TargetC', inputs=['x.txt', 'y.txt'], outputs=['z.txt']) << """
    cat x.txt y.txt > z.txt
    """

In this workflow, ``TargetA`` and ``TargetB`` each produce a file. ``TargetC``
declares that it needs two files as inputs. Since the file names match the
file names produced by `TargetA`` and ``TargetB``, ``TargetC`` depends on these
two targets.

Let's try to run this workflow:

.. code-block:: console

    $ gwf -v info run
    Scheduling target TargetC.
    Scheduling dependency TargetA of TargetC.
    Submitting target TargetA.
    Scheduling dependency TargetB of TargetC.
    Submitting target TargetB.
    Submitting target TargetC.

Notice that *gwf* first attempts to submit ``TargetC``. However, because of the
file dependencies it first schedules each dependency and submits those to the
backend. It then submits ``TargetC`` and makes sure that it will only be run
when both ``TargetA`` and ``TargetB`` has been run. If we decided that we needed
to re-run ``TargetC``, but not ``TargetA`` and ``TargetB``, we could just delete
``z.txt`` and run ``gwf run`` again. *gwf* will automatically figure out that it
only needs to run ``TargetC`` again and submit it to the backend.

What happens if we do something nonsensical like declaring a cyclic dependency?
Let's try::

    from gwf import Workflow

    gwf = Workflow()

    gwf.target('TargetA', inputs=['x.txt'], outputs=['x.txt']) << """
    echo "this is x" > x.txt
    """

Run this workflow. You should see the following:

.. code-block:: console

    error: Target TargetA depends on itself.

Observing Target Execution
--------------------------

As workflows get larger they make take a very long time to run. With *gwf* it's
easy to see how many targets have been completed, how many failed and how many
are still running using the ``gwf status`` command. We'll modify the workflow
from earlier to fake that each target takes some time to run::

    from gwf import Workflow

    gwf = Workflow()

    gwf.target('TargetA', inputs=[], outputs=['x.txt']) << """
    sleep 20 && echo "this is x" > x.txt
    """

    gwf.target('TargetB', inputs=[], outputs=['y.txt']) << """
    sleep 30 && echo "this is y" > y.txt
    """

    gwf.target('TargetC', inputs=['x.txt', 'y.txt'], outputs=['z.txt']) << """
    sleep 10 && cat x.txt y.txt > z.txt
    """

Now run ``gwf status`` (Remember to remove ``x.txt``, ``y.txt`` and ``z.txt``,
otherwise *gwf* will not submit the targets again). You should see something like this,
but with pretty colors.

.. code-block:: console

    TargetC [..............................................................] 0/0/0/3

By default, a single line is shown for each *endpoint* target. An endpoint is a target
which no other targets depends on. We can therefore see it as a kind of final target.
The dots mean that all targets in the workflow should run, but have not been submitted.
Let's try to run the workflow and see what happens.

.. code-block:: console

    $ gwf run
    $ gwf status
    TargetC [RRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS] 0/1/2/0

The ``R`` shows that one third of the targets are running (since I'm only running with one worker,
only one target can run at a time) and the other two thirds have been submitted. Running the
status command again after some time should show something like this.

.. code-block:: console

    TargetC [CCCCCCCCCCCCCCCCCCCCCRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSS] 1/1/1/0

Now the target that was running before has completed, and another target is now running,
while the final target is still just submitted. After some time, run the status command again.
All targets should now have completed, so we see this.

.. code-block:: console

    TargetC [CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC] 3/0/0/0

As you may have noticed, the numbers to the right show the number of targets that are
in a specific state in the order: completed, running, submitted, should run.

.. _function_templates:

Reusable Targets with Templates
-------------------------------

Often you will want to reuse a target definition for a lot of different files. For example,
you may have two files with reads that you need to map to a reference genome. The
mapping is the same the two files, so it would be annoying to repeat it in the workflow
specification.

Instead, *gwf* allows us to define a template which
can be used to generate one or more targets easily. In general, a template is just a function
which returns a tuple containing four values:

1. The *inputs* files,
2. The *outputs* files,
3. a dictionary with options for the target that is to be generated, for example how many
   cores the template needs and which files it depends on,
4. a string which contains the specification of the target that is to be generated.

Templates are great because they allow you to reuse functionality and encapsulate target
creation logic. Let's walk through the example above.

.. note::

    Code and data files for this example is available
    `here <https://github.com/mailund/gwf/blob/master/examples/readmapping/>`_.
    To get started, follow these steps:

    #. Change your working directory to the readmapping directory.
    #. Run ``conda env create`` to create a new environment called `readmapping`. This will
       install all required packages, including `gwf` itself, `samtools` and `bwa`.
    #. Activate the environment with ``source activate readmapping``.
    #. Open another terminal and navigate to the same directory.
    #. Activate the environment in this terminal too, using the same command as above.
    #. Start a pool with two workers with ``gwf workers -n 2``.
    #. Jump back to the first terminal. Configure *gwf* to use the local backend for this
       project using ``gwf config backend local``.
    #. You should now be able to run ``gwf status`` and all of the other *gwf* commands
       used in this tutorial.

Our reference genome is stored in ``ponAbe2.fa.gz``, so we'll need to unzip it first.
Let's write a template that unpacks files.

.. literalinclude:: ../examples/readmapping/workflow.py
   :pyobject: unzip

This is just a normal Python function that returns a tuple. The function takes two
arguments, the name of the input file and the name of the output file. In the function
we define the inputs and outputs files, a dictionary that defines the options of the targets created with this
template, and a string describing the action of the template.

We can now create a concrete target using this template::

  gwf.target_from_template('UnzipGenome',
                           unzip(inputfile='ponAbe2.fa.gz',
                                 outputfile='ponAbe2.fa'))

You could run the workflow now. The ``UnzipGenome`` target would be scheduled and submitted,
and after a few seconds you should have a ``ponAbe2.fa`` file in the project directory.

Let's now define another template for indexing a genome.

.. literalinclude:: ../examples/readmapping/workflow.py
   :pyobject: bwa_index

This template looks more complicated, but really it's the same thing as before. We define
the inputs and outputs, a dictionary with options and a string with the command that will be executed.

Let's use this template to create a target for indexing the reference genome::

  gwf.target_from_template('IndexGenome',
                           bwa_index(ref_genome='ponAbe2'))

Finally, we'll create a template for actually mapping the reads to the reference.

.. literalinclude:: ../examples/readmapping/workflow.py
   :pyobject: bwa_map

This is much the same as the previous template. Here's how we're going to use it::


  gwf.target_from_template('MapReads',
                           bwa_map(ref_genome='ponAbe2',
                                   r1='Masala_R1.fastq.gz',
                                   r2='Masala_R2.fastq.gz',
                                   bamfile='Masala.bam'))


As you can see, templates are just normal Python functions and thus they can be inspected
and manipulated in much the same way. Also, templates can be put into modules and imported
into your workflow files to facilitate reuse. It's all up to you!

Viewing Logs
------------

We may be curious about what the ``MapReads`` target wrote to the console when the target
ran, to see if there were any warnings. If a target failed, it's also valuable to see
it's output to diagnose the problem. Luckily, *gwf* makes this very easy.

.. code-block:: console

    $ gwf logs MapReads

When you run this command you'll see nothing. This is because the ``gwf logs`` command by
default only shows things written to stdout by the target, and not stderr, and apparently
nothing was written to stdout in this target. Let's try to take a look at stderr instead
by applying the ``--stderr`` flag.

.. code-block:: console

    $ gwf logs --stderr MapReads
    [M::bwa_idx_load_from_disk] read 0 ALT contigs
    [M::process] read 15000 sequences (1500000 bp)...
    [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (1, 65, 1, 0)
    [M::mem_pestat] skip orientation FF as there are not enough pairs
    [M::mem_pestat] analyzing insert size distribution for orientation FR...
    [M::mem_pestat] (25, 50, 75) percentile: (313, 369, 429)
    [M::mem_pestat] low and high boundaries for computing mean and std.dev: (81, 661)
    [M::mem_pestat] mean and std.dev: (372.88, 86.21)
    [M::mem_pestat] low and high boundaries for proper pairs: (1, 777)
    [M::mem_pestat] skip orientation RF as there are not enough pairs
    [M::mem_pestat] skip orientation RR as there are not enough pairs
    [M::mem_process_seqs] Processed 15000 reads in 1.945 CPU sec, 0.678 real sec
    [main] Version: 0.7.15-r1140
    [main] CMD: bwa mem -t 16 ponAbe2 Masala_R1.fastq.gz Masala_R2.fastq.gz
    [main] Real time: 0.877 sec; CPU: 2.036 sec

We can do this for any target in our workflow. The logs shown are always the most recent
ones since *gwf* does not archive logs from old runs of targets.

Cleaning Up
-----------

Now that we have run our workflow we may wish to remove intermediate files to save disk
space. In *gwf* we can use the ``gwf clean`` command for this. If we wanted to remove all
files produced by all targets in the workflow, we could just run:

.. code-block:: console

    $ gwf clean

However, this would also remove the file with the mapped genome, which we probably want
to save. This file was produced by an endpoint target (a target which no other
target depends on), so we can tell *gwf* to only remove files that were not produced
by an endpoint:

.. code-block:: console

    $ gwf clean --not-endpoints
    INFO    |  Cleaning up: IndexGenome.
    INFO    |  Cleaning up: UnzipGenome.


A Note About Reproducibility
----------------------------

Reproducibility is an important part of research and since *gwf* workflows describe every
step of your computation, how the steps are connected, and the files produced in each step,
it's a valuable tool in making your workflows reproducible. In combination with the
``conda`` package manager and the concept of environments, you can build completely
reproducible workflows in a declarative, flexible fashion.

Consider the read mapping example used above. Since we included a specification of the
complete environment through a ``environment.yml`` file, which even included samtools,
bwa and *gwf* itself, we were able to easily create a working environment with exactly
the right software versions used for our workflow. The whole workflow could also easily
be copied to a cluster and run through e.g. the Slurm backend, since we can exactly
reproduce the environment used locally.


.. _Anaconda: https://www.continuum.io/downloads
