.. _writing_backends:

================
Writing Backends
================

Backends in *gwf* are the interface between *gwf* and whatever can be used to
execute a target. For example, the Slurm backend included with *gwf* submits
targets to the `Slurm Workload Manager`_.

To get started we must declare define a class that inherits from
:class:`~gwf.backends.Backend`::

    # mybackend/mybackend.py
    from gwf.backends import Backend


    class MyBackend(Backend):
        pass

Registering a Backend
=====================

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

Backends must implement a set of methods that *gwf* uses to submit to the backend
and query the backend for the status of targets.

Option Defaults
===============

The backend should define :attr:`option_defaults` as an attribute. The value must
be a dictionary mapping option names to defaults, e.g.::

    # mybackend/mybackend.py
    from gwf.backends import Backend

    class MyBackend(Backend):
        option_defaults = {
            'cores': 1,
        }

Internally, *gwf* uses this dictionary to check whether targets contain options
not supported by the backend and warn the user if this is the case. Thus, *all*
options supported by the backend must be specified in this dictionary.

Targets in a workflow can now declare the number of cores they wish to allocate
and we can use this information in :func:`~gwf.backends.Backend.submit` to allocate the
given number of cores for the target to be submitted. If the user doesn't specify
the `cores` option it will default to 1.

If want to specify support for an option, but there is no sensible default value
(e.g. in the case of a username or e-mail address), use `None` as the default value.

Implementing the Backend Interface
==================================

Our backend still doesn't really do anything. We've only told *gwf* that our backend
exists (by its entrypoint) and which options are supported. To get a backend to
actually work we must implement three methods: :func:`~gwf.backends.Backend.submit`,
:func:`~gwf.backends.Backend.cancel` and :func:`~gwf.backends.Backend.status`. If needed,
one may also implement the :func:`~gwf.backends.Backend.close` method, which will be
called when the backend is no longer needed (right before *gwf* exits).

All methods must return immediately, that is, calling :func:`~gwf.backends.Backend.submit`
should submit the target for execution in some other process, but not run the target
itself. For example, the local backend connects to a set of workers running in
a different process, submits jobs to these workers and returns immedicately.


Storing Log Files
=================

Backends can store log files in different ways. For example, the Slurm backend stores log
files as files on disk, while other backends may wish to store log files in an S3 bucket
or in a database.

To allow for all of these scenarios, *gwf* has the concept of a *log manager*. The log
manager interface only assumes that log files can be written and accessed through file-like
objects. Log managers should inherit from :class:`~gwf.backends.base.LogManager`.

.. code-block:: python

    from gwf.backends.base import LogManager

    class MyLogManager(LogManager):

        def open_stdout(self, target, mode='r'):
            pass

        def open_stderr(self, target, mode='r'):
            pass

Each method must return a file-like object providing access to the log data for `target`.
Log managers can also provide other methods. For example, the
:class:`~gwf.backends.base.FileLogManager` provides methods for retrieving the paths of the
log files.

Backends should set the :attr:`log_manager` attribute on the class to an instance of a the
log manager to be used. The log manager must be set as a class attribute to allow access to
log files without initializing the backend, which may be slow.

At the moment we provide two log managers:

* :class:`~gwf.backends.base.FileLogManager` (default)
* :class:`~gwf.backends.base.MemoryLogManager`


Handling Configuration
======================

We can allow the user to configure aspects of the backend by using the central
configuration object.

.. code-block:: python

    from gwf.conf import config

    key1 = config.get('yourbackend.key1', 'default1')
    key2 = config.get('yourbackend.key2', 'default2')

Backends should provide reasonable defaults, as shown above.
The user can set configuration keys using the builtin ``config`` command::

    $ gwf config set yourbackend.key1 value1
    $ gwf config set yourbackend.key2 value2


.. _Slurm Workload Manager: http://slurm.schedmd.com/
