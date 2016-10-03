.. _writing_plugins:

Writing Plugins
===============

Plugins can manipulate and control the entire workflow and its execution. Thus,
quite complex things can be implemented with plugins. Actually, much of the
main functionality of *gwf* is implemented through plugins! In this tutorial
we will implement a plugin that prints the names of all targets in the workflow.

A plugin is a Python package and thus you should be (at least slightly)
familiar with the principles of packaging in Python. If you're not, don't worry.
This guide will teach you the basics.

Plugins for *gwf* must inherit from :class:`gwf.plugins.Plugin` and
implement the required methods and specify the :attr:`name` attribute. Thus,
a minimal (and completely useless) plugin looks like this::

    # myplugin/myplugin.py
    from gwf.plugins import Plugin

    class MyPlugin(Plugin):
        name = 'myplugin'

This plugin does nothing, but let's try to hook it up to *gwf*. Plugins are
registered with *gwf* through entry points specified in `setup.py`, so we have
to write a `setup.py` file for our plugin::

    # myplugin/setup.py
    from setuptools import setup

    setup(
        name="myplugin",
        version="0.0.1",
        py_modules=['myplugin'],
        install_requires=[
          'gwf>=1.0',
        ],
        entry_points={
            'gwf.plugins': [
                'myplugin = myplugin:MyPlugin',
            ]
        },
    )

Here we register our plugin under the `gwf.plugins` entrypoint. This means that
*gwf* will be able to automatically locate and load our plugin (*gwf* looks
for all entry points registered under the `gwf.plugins` identifier). Since the
plugin is written for *gwf*, we also specify that our package requires at least
version 1.0 of *gwf*.

We can now install our plugin by running `python setup.py develop`. We install
the plugin in development mode so that we can quickly test changes.

When the installation has completed we can check whether *gwf* can find and
load the plugin by running `gwf -b testing -v debug` (you must do this in a
directory containing a `workflow.py` file containing a valid workflow).

You should see something like this::

    DEBUG |  Platform: Darwin-16.0.0-x86_64-i386-64bit.
    DEBUG |  GWF version: 1.0.0.
    DEBUG |  Python version: 3.5.2.
    DEBUG |  Node: Dans-Mac-mini.local.
    DEBUG |  Python path: ['/Users/das/Code/gwf/examples/minimal-workflow', ...]
    DEBUG |  Loaded plugins: logs, run, config, myplugin.
    DEBUG |  Loaded backends: testing, slurm.
    DEBUG |  Setting active backend: testing.
    DEBUG |  Loading workflow from: examples/minimal-workflow/workflow.py.
    ...

Some of the output has been cut away, but you should be able to locate
`myplugin` in the list of loaded plugins.

.. note::
  The observant reader may have noticed that there are several plugins included
  with *gwf*, namely the `logs`, `run` and `config` plugins that provide the
  `logs`, `run` and `config` commands, respectively. In short, we eat our own
  dog food!

So now our plugin is loaded by *gwf*, but it doesn't do anything. Let's set up
a subcommand that the user can run to get a list of all targets in the
current workflow. Plugins have a special hook to set up the argument parsing:

.. automethod::
  gwf.plugins.Plugin.setup_argument_parser
  :noindex:

We can override this method to add our subcommand as follows::

    # myplugin/myplugin.py
    from gwf.plugins import Plugin

    class MyPlugin(Plugin):
        name = 'myplugin'

        def setup_argument_parser(self, parser, subparsers):
            subparser = self.setup_subparser(
                subparsers,
                'print-targets',
                'A command for printing the name of all targets.',
                self.on_run,
            )

        def on_run(self):
            print('hello!')

Here we use the helper method :func:`setup_subparser` to
hook up the ``on_run`` method to the subcommand. Running ``gwf -h`` should now
list the ``print-targets`` subcommand and running ``gwf print-targets`` should
print ``hello!`` to the screen.

Now we just need to get a list of all targets in the workflow and print their
names. Plugins are by default configured with three attributes: :attr:`backend`,
:attr:`config`, and :attr:`workflow`. The :attr:`workflow` attribute contains a
:class:`~gwf.PreparedWorkflow` which has a :attr:`targets` attribute. We can use
this to look up all targets in the workflow::

    # myplugin/myplugin.py
    from gwf.plugins import Plugin

    class MyPlugin(Plugin):
        name = 'myplugin'

        def setup_argument_parser(self, parser, subparsers):
            subparser = self.setup_subparser(
                subparsers,
                'print-targets',
                'A command for printing the name of all targets.',
                self.on_run,
            )

        def on_run(self):
            for target_name in self.workflow.targets.keys():
                print(target_name)


Voilà! We now have a functioning plugin that prints a list of all targets in
the workflow.

Adding and Using Command Options
--------------------------------


Logging from a Plugin
---------------------

Write something about logging...


Examples
--------

The plugins included with *gwf* are good examples of the kind of functionality
that can be implemented through plugins. You can find them
`here <https://github.com/mailund/gwf/tree/master/gwf/plugins>`_.

Further Reading
---------------

Link to Python Packaging guide.
