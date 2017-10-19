===============
Large Workflows
===============

While *gwf* can handle quite large workflows without any problems, there are some things that may cause significant pain
when working with very, very large workflows, especially when the workflows has many (> 50000) targets producing many
files. However, the problems depend hugely on your filesystem since most scalability problems are caused by the time it
takes *gwf* to access the filesystem when scheduling targets.

In this section we will show a few tricks for handling very large workflows.

I have to run the same pipeline for *a lot* of files and running ``gwf status`` is very slow.
---------------------------------------------------------------------------------------------

In this case *gwf* is probably slow because computing the dependency graph for your entire workflow takes a while and
because *gwf* needs to access the filesystem for each input and output file in the workflow to check if any targets
should be re-run.

One solution to this problem is to dynamically generate individual workflows for each input file, as shown here:

.. code-block:: python

    from glob import glob
    from gwf import Workflow

    data_files = ['Sample1', 'Sample2', 'Sample3']
    for input_file in data_files:
        workflow_name = 'Analyse.{}'.format(input_file)

        wf = Workflow(name=workflow_name)
        wf.target('{}.Filter'.format(input_file), inputs=[input_file], outputs=[...]) << """..."""
        wf.target('{}.ComputeSummaries'.format(input_file), ...) << """..."""

        globals()[workflow_name] = wf

You can now run the workflow for a single sample by specifying the name of the workflow:

.. code-block:: console

    $ gwf run -f workflow.py:Analyse.Sample1 run

This will only run the targets associated with `Sample1`. While this means that running *all* workflows in one go
involves a bit more work, it also means that *gwf* will only have to compute the dependency graph and check timestamps
for the targets associated with the selected sample.
