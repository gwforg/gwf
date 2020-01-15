==========
Change Log
==========


Version 1.7.1
=============

Changed
-------

* The ``clean`` command now asks for confirmation when no targets are
  specified. This can be avoided by using the ``--force`` flag.
* The Slurm backend will now get all jobs from ``squeue`` and not only jobs
  belonging to the current user. This used to be the default behaviour, but was
  unintentionally changed when refactoring the Slurm backend.

Fixed
-----

* The output of the ``cancel`` command is now correct.


Version 1.7.0
=============

This release is the first version of *gwf* to support Python 3.8. We're now
testing on nightly version of Python so that future Python releases should be
supported as soon as they're available.

Added
-----

* Automatically checks for updates on a regular basis. A warning will be shown
  when a new version is available. This feature can be disabled with
  ``gwf config set check_updates no``.
* The Slurm backend will now retry failed operations, e.g. if it fails to get
  the state of the queue from Slurm or if submitting a job failed.
* Support for Python 3.8.

Changed
-------

* When running `gwf cancel` without specifying targets, only
  submitted/running targets will be cancelled.

Fixed
-----

* Fixed a deadlock in the local backend which caused the workers to freeze.
* Links to examples and ReadTheDocs (we're now hosting the docs ourselves).


Version 1.6.0
=============

This release contains a few large features such as named inputs and outputs
and the introduction of the :class:`Workflow.map` method for easily generating
multiple targets from a template and improving readability.

There are also several minor improvements that contribute to the overall
user experience ranging from speed improvements for the ``logs`` command to
improved debugging output.

The documentation has also been restructured and improved.

Added
-----

* Named inputs and outputs. The ``inputs`` and ``outputs`` arguments to
  ``Workflow.target`` can now be either a string, list or dictionary. See the
  documentation for more details.
* Tutorial now explains what happens if a target fails.
* Documentation now has an official list of *gwf* plugins.
* The ``status`` command now has a ``--summary`` option the summarizes the
  status of  an entire workflow.
* All input and output paths are now checked for non-printable characters such
  as newlines and tabs. This can cause problems that are very hard to find and
  fix (e.g. *gwf* reporting that a file is missing even though it seems to be
  there). Paths containing such characters now result in an error.

Changed
-------

* The ``logs`` command is now much faster since it no longer builds the entire
  graph.
* The target in the ``status`` output is now sorted in creation order, instead
  of alphabetically.
* Cleaner output formatting, especially when running with ``-v debug``.
* Improved log messages from scheduler when running with ``-v debug``. The
  messages are now more specific and helpful.
* Documentation has been restructured to be more readable and have less
  redundant information.

Fixed
-----

* Crash when running ``gwf init`` without an existing configuration file.
* The ``--force`` flag for the ``cancel`` command now actually forces cancellation.
* We now respect the ``--no-color`` flag completely and implement the ``NO_COLOR``
  environment variable standard described `here <https://no-color.org/>`_.


Version 1.5.1
=============

Fixed
-----

* Crash when Slurm returns unknown job state (#244).

Version 1.5.0
=============

Added
-----

* Users can now run ``gwf init`` to bootstrap a new *gwf* project (c78193).
* Add option to protect output files in a target from being removed when
  ``gwf clean`` is being run (2f51ed).

Fixed
-----

* Ensure job script end with a newline (#239).
* Ignore missing log files when cleaning on run (#237).

Version 1.4.0
=============

Added
-----

* Backend for Sun Grid Engine (SGE). The backend does not support all target
  options supported by the Slurm backend, so workflows can not necessarily
  run with the SGE backend without changes. See the documentation for a list
  of supported options.

Version 1.3.2
=============

Fixed
-----

* Made the ``touch`` command faster.

Version 1.3.1
=============

Added
-----

* The ``gwf status`` command now accepts multiple ``-s/--status`` flags and will show
  targets matching any of the given states. E.g. ``gwf status -s completed -s running``
  will show all completed and running targets.
* A new command ``gwf touch`` has been introduced. The command touches all files in
  the workflow in order, creating missing files and updating timestamps, such that
  *gwf* thinks that the workflow has been run.
* When specifying the workflow attribute in the workflow path, e.g.
  ``gwf -f workflow.py:foo``, the filename part can now be left out and will default
  to `workflow.py`. For example, ``gwf -f :foo`` will access the ``foo`` workflow
  object in `workflow.py`.
* Documentation describing advanced patterns for *gwf* workflows.


Version 1.3.0
=============

This release contains a bunch of new features and plenty of bug fixes. Most
noteworthy is the removal of the progress bars in the status command. The status
bars were often confusing and didn't communicate much more than a simple
"percentage completion". The status command now outputs a table with target
name, target status, and percentage completion (see the tutorial for examples).
Additionally, the status command now shows all targets by default (not only
endpoints). For users who wish to only see endpoints, there's now a
``--endpoints`` flag.

We aim to make *gwf* a good cluster citizen. Thus, logs from targets that no
no longer exist in the workflow will now be removed when running ``gwf run``.
This ensures that *gwf* doesn't unnecessarily accumulate logs over time.

Fixed
-----

* Add missing import to documentation for function templates (4eddcac).
* Remove reference to ``--not-endpoints`` flag (d7ed251).
* Remove broken badges in README (e352f09).
* Remove pre-1.0 upgrade documentation (bfa03da6).
* Fixed bug in scheduler that caused an exception when a target's input file did
  not exist, but the output file did (reported by Jonas Berglund) (92301ef3).

Changed
-------

* Dots have been removed from logging output to make copy-pasting target names
  easier (f33f7195).
* Now uses pipenv to fix development environment.
* Improved coloring of logging output when running with ``-v debug`` (ab4ac7e3).
* Remove status bars in ``gwf status`` command (47cb7b50).

Added
-----

* Added undocumented API which allows core and plugins to register validation
  functions for configuration keys. This fixes issues like #226 (c8c57d7c7).
* The ``gwf clean`` command now shows how much data will be removed (d81f143f1).
* Remove log files for targets that are no longer defined in the workflow
  (beb912bd).
* Note in tutorial on how to terminate the local workers (a long with other
  updates to the tutorial) (34421498).

Version 1.2.1
=============

Fixed
-----

* Bug when returning an ``AnonymousTarget`` from a template function without
  specifying the *working_dir* in the constructor (#212). Thanks to Steffen
  Møller-Larsen for reporting this.

Version 1.2
===========

Fixed
-----

* Bug when using ``--format table`` and no targets were found (#203).
* Bug when cancelling a target running on the Slurm backend (#199).
* Link to documentation in error message when unable to connect to local
  workers.
* Fixed bug in the *FileLogManager* where the wrong exception was raised when no
  log was found.

Changed
-------

* Moved checking of file timestamps to the scheduler. This means that creating a
  ``Graph`` object will never touch the file system, and thus won't raise an
  exception if a target depends on a file that doesn't exist and that's not
  provided a target. Instead, unresolved paths are added to
  ``Graph.unresolved``. They will then be checked by the scheduler (if
  necessary). For end users, this means that many commands have become
  substantially faster.

Added
-----

* Added ``AnonymousTarget`` which represents an unnamed target. ``Target`` now
  inherits from this class and templates may now return an ``AnonymousTarget``
  instead of a tuple.
* Added *backend.slurm.log_mode* option, see the documentation for the Slurm
  backend for usage (#202).

Version 1.1
===========

Fixed
-----

* Very slow scheduling when using dry run with unsubmitted targets (#184, 93e71a).
* Fixed cancellation with the Slurm backend (#183, 29445f).
* Fixed wildcard filtering of targets (#185, 036e3d).

Changed
-------

* Move file cache construction out of ``Graph`` (#186, 93e71a). This change is
  invisible to end-users, but speeds up the ``logs``, ``cancel``, ``info``,
  ``logs`` and ``workers`` commands.
* Replaced ``--not-endpoints`` flag in ``clean`` command with ``--all`` flag.
* Made filtering more intuitive in all commands.
* The ``info`` command now outputs JSON instead of invalid YAML.
* The ``info`` command outputs information for all targets in the workflow by
  default.
* Backends must now specify a ``log_manager`` class attribute specifying which
  log manager to use for accessing target log files.
* Backends should now be used as context managers to make sure that
  ``Backend.close()`` is called when the backend is no longer needed, as it is
  no longer called automatically on exit.

Added
------

* Added filtering of targets by name in the ``info`` command.
* Added API documentation for the ``gwf.filtering`` module.
* Added ``gwf.core.graph_from_path()`` and ``gwf.core.graph_from_config()``.
* Added ``gwf.backends.list_backends()``, ``gwf.backends.backend_from_name()``
  and ``gwf.backends.backend_from_config()``.
* Added ``SlurmBackend.get_job_id()`` and ``SlurmBackend.forget_job()`` to
  ``SlurmBackend`` to make it easier for plugins to integrate with Slurm.
* Documentation for log managers.
* Documentation on how to handle large workflows.


Version 1.0
===========

First stable release of *gwf*! We strongly encourage users of pre-1.0 users to
read the tutorial, since quite a lot of things have changed. We also recommend
reading the guide for converting pre-1.0 workflows to version 1.0. However,
users attempting to do this should be aware that the the template mechanism in
1.0 is slightly different and thus requires rewriting template functions.

Fixed
-----

* Fixed a bug which caused *gwf* to fail when cancelling jobs when using the
  Slurm backend (8c1717).

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

* Fixed a subtle bug in scheduling which caused problems when resubmitting a
  workflow where some targets were already running (a5d884).
* Fixed a bug in the ``SlurmBackend`` which caused *gwf* to crash if the Slurm
  queue contained a job with many dependencies (eb4446).
* Added back the `-e` flag in the ``logs`` command.


Version 1.0b9
=============

Fixed
-----

* Fixed a bug in the ``SlurmBackend`` which caused running targets as unknown
  (33a6bd).

Changed
-------

* The Slurm backend's database of tracked jobs is now cleaned on initialization
  to keep it from growing indefinitely (bd3f95).

Version 1.0b8
=============

Fixed
-----

* Fixed a bug which caused the *gwf logs* command to always show stderr
  (01b267).

* Fixed a bug which caused dependencies to be set incorrectly when two targets
  depended on the same target (4d9e07).

Changed
-------

* Improved error message when trying to create a target from an invalid template
  (d27d1f).

* Improved error message when assigning a non-string spec to a target (2aca0a).

* `gwf logs` command now outputs logs via a pager when the system supports it,
  unless `--no-pager` is used (01b267).

Added
-----

* Added more tests to cover scenarios with included workflows when building the
  workflow graph (86a68d0).

* Added a bunch of documentation (69e136, 51a0e7, 942b05).

Version 1.0b7
=============

Fixed
-----

* Fixed bug in scheduling which was actually the cause of the incorrect
  scheduling that was "fixed" in 1.0b6. Also added documentation for
  ``gwf.core.schedule`` (7c47cb).

Changed
-------

* Updated documentation in a bunch of places, mostly styling.

Version 1.0b6
=============

Fixed
-----

* A bug in ``SlurmBackend`` which caused dependencies between targets to not be
  set correctly (6b71d2).

Changed
-------

* More improvements to and clean up of build process.
* Updated some examples in the tutorial with current output from *gwf* (42c5da).
* Logging output is now more consistent (b95af04).

Added
-----

* Documentation for maintainers on how to merge in contributions and rolling a
  new release (fe1ee3).

Version 1.0b5
=============

Fixed
-----

* Unset option passed to backend causes error (#166, dcff44).
* Set import path to allow import of module in workflow file (64841c).

Changed
-------

* Vastly improved build and deploy process. We're now actually building and
  testing with conda.
