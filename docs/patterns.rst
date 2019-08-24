.. _patterns:

========
Patterns
========

This guide takes you through some advanced features and patterns that can
be utilized in *gwf*. Remember that *gwf* is just a way of generating
workflows using the Python programming language and thus many of these
patterns simply use plain Python code to abstract and automate certain
things.

Iterating Over a Parameter Space
--------------------------------

Say that you have a workflow that runs a program with many different
combinations of parameters, e.g. the parameters *xs*, *ys*, and *zs*. Each
parameter can take multiple values:

.. code-block::: python

    xs = [0, 1, 2, 4, 5]
    ys = ['cold', 'warm']
    zs = [0.1, 0.2, 0.3, 0.4, 0.5]

We now want to run out program `simulate` with all possible combinations
of these parameters. To do this, we'll use the Python function
:func:`itertools.product()` to create an iterator over all combinations
of the parameters:

.. code-block:: python

    import itertools

    parameter_space = itertools.product(xs, ys, zs)

We can then iterate over the parameter space:

.. code-block:: python

    gwf = Workflow()

    for x, y, z in parameter_space:
        gwf.target(
            name='sim_{}_{}_{}'.format(x, y, z),
            inputs=['input.txt'],
            outputs=['output_{}_{}_{}.txt'.format(x, y, z)],
        ) << """
        ./simulate {} {} {}
        """.format(x, y, z)


Dynamically Generating a Workflow
---------------------------------

We can make our workflows more reusable by generating them dynamically. For
example, we may wish to make it easy for others to change the inputs to our
workflow or let users specify a different output directory. When generating
workflows dynamically you can essentially parameterize the workflow in any
way you want. In combination with inclusion of workflows into other
workflows, this allows for extremely powerful composition.

To dynamically generate a workflow, we simply create a function which
builds the workflow and returns it:

.. code-block:: python

    import os.path.join
    from gwf import Workflow

    def my_fancy_workflow(output_dir='outputs/'):
        # Create an empty workflow object.
        w = Workflow()

        # Add targets to the workflow object, respecting the value of `output_dir`.
        foo_output = os.path.join(output_dir, 'output1.txt')
        w.target(
            name='Foo',
            inputs=['input.txt'],
            outputs=[foo_output],
        ) << """
        ./run_foo > {}
        """.format(foo_output)

        bar_output = os.path.join(output_dir, 'output2.txt')
        w.target(
            name='Bar',
            inputs=[foo_output],
            outputs=[bar_output]
        )

        # Now return the workflow.
        return w


You can put this function in file next to your workflow, or any other place from
which you can import the function. In this case, let's put the file next to
``workflow.py`` in a file called ``fancy.py``.

In ``workflow.py`` we can then use the workflow as follows:

.. code-block:: python

    from fancy import my_fancy_workflow

    gwf = my_fancy_workflow()

We can now run the workflow as usual:

.. code-block:: shell

    $ gwf run

However, we can now easily change the output directory:

.. code-block:: python

    from fancy import my_fancy_workflow

    gwf = my_fancy_workflow(output_dir='new_outputs/')

Parameterizing the workflow can also let the user choose to deactivate parts of
the workflow. For example, imagine that ``Bar`` generates summary files that may
now always be needed. In this case, we can let the user choose to leave it out:

.. code-block:: python

    import os.path.join
    from gwf import Workflow

    def my_fancy_workflow(output_dir='outputs/', summarize=True):
        # Create an empty workflow object.
        w = Workflow()

        # Add targets to the workflow object, respecting the value of `output_dir`.
        foo_output = os.path.join(output_dir, 'output1.txt')
        w.target(
            name='Foo',
            inputs=['input.txt'],
            outputs=[foo_output],
        ) << """
        ./run_foo > {}
        """.format(foo_output)

        # Only create target `Bar` if we want to summarize the data.
        if summarize:
            bar_output = os.path.join(output_dir, 'output2.txt')
            w.target(
                name='Bar',
                inputs=[foo_output],
                outputs=[bar_output]
            )

        # Now return the workflow.
        return w

In ``workflow.py`` we can then use the workflow as follows:

.. code-block:: python

    from fancy import my_fancy_workflow

    gwf = my_fancy_workflow(summarize=False)


External Configuration of Workflows
-----------------------------------

In the previous section we saw how we can parameterize workflows. However, in some
cases we may want to let the user of our workflow specify the parameters without
touching any Python code at all. That is, we want an external configuration file.

The configuration format could be anything, but in this example we'll use a JSON
as the configuration format. First, this is what our configuration file is going
to look like:

.. code-block:: json

    {
        "output_dir": "some_output_directory/",
        "summarize": true
    }

We put this file next to ``workflow.py``, e.g. as ``config.json``. We can now read
the configuration using the Python ``json`` module in ``workflow.py``:

.. code-block:: python

    import json
    from fancy import my_fancy_workflow

    config = json.load(open('config.json'))

    gwf = my_fancy_workflow(
        output_dir=config['output_dir'],
        summarize=config['summarize'],
    )

We can now change the values in ``config.json`` and run the workflow as usual.
