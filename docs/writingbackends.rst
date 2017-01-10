.. _writing_backends:

Writing Backends
================

Backends in *gwf* are the interface between *gwf* and whatever can be used to
execute a target.

To implement a backend you should first read
:ref:`writing_plugins` since backends are implemented and registered
much the same way as plugins.

To get started we must declare define a class that inherits from
:class:`gwf.backends.Backend`::

    # mybackend/mybackend.py
    from gwf.backends import Backend

    class MyBackend(Backend):
        name = 'mybackend'

Registering a Backend
---------------------

Backends must be registered under the ``gwf.backends`` entrypoint as shown
here::

    # mybackend/setup.py
    from setuptools import setup

    setup(
        name="mybackend",
        version="0.0.1",
        py_modules=['mybackend'],
        install_requires=[
          'gwf>=1.0',
        ],
        entry_points={
            'gwf.backends': [
                'mybackend = mybackend:MyBackend',
            ]
        },
    )

A backend can manipulate *gwf* in much the same way a plugin can. For example,
the backend could add an argument to the command line interface, add
subcommands etc. However, backends must implement a set of methods that *gwf*
uses to submit and query targets.

Specifying Target Options
-------------------------

The backend must define the two class variables
:attr:`supported_options` and :attr:`option_defaults`. The first attribute
declares which options the backend supports. For example, we may want the user
to be able specify the number of cores to be allocated for a given target::

    # mybackend/mybackend.py
    from gwf.backends import Backend

    class MyBackend(Backend):
        name = 'mybackend'
        supported_options = ['cores']

Targets in a workflow can now declare the number of cores they wish to allocate
and we can use this information in :func:`Backend.submit` to allocate the
given number of cores for the target to be submitted. However, the user should
not be forced to specify the number of cores for every target in the workflow,
so we should supply a sensible default. This is accomplished through
:attr:`option_defaults`::

    # mybackend/mybackend.py
    from gwf.backends import Backend

    class MyBackend(Backend):
        name = 'mybackend'
        supported_options = ['cores']
        option_defaults = {
            'cores': 1,
        }

Now all targets will by default have 1 core allocated unless the user has
specified something else for the target.



Reference
---------

.. autoclass:: gwf.backends.Backend
   :members: supported_options, option_defaults, submit, cancel, submitted,
      running, logs
   :member-order: bysource

.. _Slurm Workload Manager: http://slurm.schedmd.com/
