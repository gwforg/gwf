==========
Change Log
==========


Version 1.2
===========

Fixed
-----

* Bug when using ``--format table`` and no targets were found (#203).
* Bug when cancelling a target running on the Slurm backend (#199).
* Link to documentation in error message when unable to connect to local workers.

Changed
-------

* Moved checking of file timestamps to the scheduler. This means that creating a ``Graph`` object will never touch the
  file system, and thus won't raise an exception if a target depends on a file that doesn't exist and that's not
  provided a target. Instead, unresolved paths are added to ``Graph.unresolved``. They will then be checked by the
  scheduler (if necessary).


Added
-----

* Added ``AnonymousTarget`` which represents an unnamed target. ``Target`` now inherits from this class and
  templates may now return an ``AnonymousTarget`` instead of a tuple.
* Added *backend.slurm.log_mode* option, see the documentation for the Slurm backend for usage (#202).

Version 1.1
===========

Fixed
-----

* Very slow scheduling when using dry run with unsubmitted targets (#184, 93e71a).
* Fixed cancellation with the Slurm backend (#183, 29445f).
* Fixed wildcard filtering of targets (#185, 036e3d).

Changed
-------

* Move file cache construction out of ``Graph`` (#186, 93e71a). This change is invisible to end-users, but speeds up the
  ``logs``, ``cancel``, ``info``, ``logs`` and ``workers`` commands.
* Replaced ``--not-endpoints`` flag in ``clean`` command with ``--all`` flag.
* Made filtering more intuitive in all commands.
* The ``info`` command now outputs JSON instead of invalid YAML.
* The ``info`` command outputs information for all targets in the workflow by default.
* Backends must now specify a ``log_manager`` class attribute specifying which log manager to use for accessing
  target log files.
* Backends should now be used as context managers to make sure that ``Backend.close()`` is called when the backend is no
  longer needed, as it is no longer called automatically on exit.

Added
------

* Added filtering of targets by name in the ``info`` command.
* Added API documentation for the ``gwf.filtering`` module.
* Added ``gwf.core.graph_from_path()`` and ``gwf.core.graph_from_config()``.
* Added ``gwf.backends.list_backends()``, ``gwf.backends.backend_from_name()`` and
  ``gwf.backends.backend_from_config()``.
* Added ``SlurmBackend.get_job_id()`` and ``SlurmBackend.forget_job()`` to ``SlurmBackend`` to make it easier for
  plugins to integrate with Slurm.
* Documentation for log managers.
* Documentation on how to handle large workflows.


Version 1.0
===========

First stable release of *gwf*! We strongly encourage users of pre-1.0 users to read the tutorial, since quite a lot
of things have changed. We also recommend reading the guide for converting pre-1.0 workflows to version 1.0. However,
users attempting to do this should be aware that the the template mechanism in 1.0 is slightly different and thus
requires rewriting template functions.

Fixed
-----

* Fixed a bug which caused *gwf* to fail when cancelling jobs when using the Slurm backend (8c1717).

Changed
-------

* Documentation in various places, especially the core API.
* Documentation for maintainers.

Added
-----

* Topic guide covering templates (b175fe).
* Added ``info`` command (6dbdbb).


Version 1.0b10
==============

Fixed
-----

* Fixed a subtle bug in scheduling which caused problems when resubmitting a workflow where some targets were already running (a5d884).

* Fixed a bug in the ``SlurmBackend`` which caused *gwf* to crash if the Slurm queue contained a job with many dependencies (eb4446).

* Added back the `-e` flag in the ``logs`` command.


Version 1.0b9
=============

Fixed
-----

* Fixed a bug in the ``SlurmBackend`` which caused running targets as unknown (33a6bd).

Changed
-------

* The Slurm backend's database of tracked jobs is now cleaned on initialization to keep it from growing indefinitely (bd3f95).

Version 1.0b8
=============

Fixed
-----

* Fixed a bug which caused the *gwf logs* command to always show stderr (01b267).

* Fixed a bug which caused dependencies to be set incorrectly when two targets depended on the same target (4d9e07).

Changed
-------

* Improved error message when trying to create a target from an invalid template (d27d1f).

* Improved error message when assigning a non-string spec to a target (2aca0a).

* `gwf logs` command now outputs logs via a pager when the system supports it, unless `--no-pager` is used (01b267).

Added
-----

* Added more tests to cover scenarios with included workflows when building the workflow graph (86a68d0).

* Added a bunch of documentation (69e136, 51a0e7, 942b05).

Version 1.0b7
=============

Fixed
-----

* Fixed bug in scheduling which was actually the cause of the incorrect scheduling that was "fixed" in 1.0b6.
  Also added documentation for ``gwf.core.schedule`` (7c47cb).

Changed
-------

* Updated documentation in a bunch of places, mostly styling.

Version 1.0b6
=============

Fixed
-----

* A bug in ``SlurmBackend`` which caused dependencies between targets to not be set correctly (6b71d2).

Changed
-------

* More improvements to and clean up of build process.
* Updated some examples in the tutorial with current output from *gwf* (42c5da).
* Logging output is now more consistent (b95af04).

Added
-----

* Documentation for maintainers on how to merge in contributions and rolling a new release (fe1ee3).

Version 1.0b5
=============

Fixed
-----

* Unset option passed to backend causes error (#166, dcff44).
* Set import path to allow import of module in workflow file (64841c).

Changed
-------

* Vastly improved build and deploy process. We're now actually building and testing with conda.
