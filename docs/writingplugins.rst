.. _writing_plugins:

Writing Plugins
===============

Plugins for *gwf* must inherit from :class:`gwf.plugins.base.Plugin` and
implement the required methods and specify the :attr:`name` attribute. Thus,
a minimal (and completely useless) plugin looks like this::

    # myplugin/myplugin.py
    from gwf.plugins.base import Plugin

    class MyPlugin(Plugin):
        name = 'myplugin'

This plugin does nothing, but let's try to hook it up to *gwf*. Plugins are
registered with *gwf* through entrypoints specified in `setup.py`, so we have
to write a `setup.py` file for our plugin::

    # setup.py
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

Here we register our plugin under the `gwf.plugins` entrypoint. Since the
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

Some of the output has been cut away, but you should be able to locate `myplugin`
in the list of loaded plugins. The observant reader may have noticed that
there are several plugins included with *gwf*, namely the `logs`, `run` and
`config` plugins that provide the `logs`, `run` and `config` commands,
respectively.
