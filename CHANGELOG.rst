==========
Change Log
==========

Version 2.1.1
=============

Fixed
-----

* The exec component would sometimes fail on Python 3.9 because it wasn't
  catching the right exception.
* The exec component would sometimes issue an empty write to stdout/stderr
  which (at least in theory) could confuse the consumer of those streams.
* Fixed syntax errors in the ``gwf run`` and ``gwf status`` commands such
  that they now work in Python 3.8 again.


Added
-----

* We now test correctly on all supported Python versions.

Version 2.1.0
=============

This is a pretty substantial release (but entirely backwards compatible).
There's new functionality submitted by several new contributors, as well as a
new component that enables nice runtime (as in, when a target is being run by
the backend) features.

Added
-----

* This release introduces the concept of *executors* which can be used to enable
  runtime behavior for targets. This means that you can now use executors to
  run your target inside Conda, Pixi, Apptainer, and Singularity environments:

  .. code-block:: python

    from gwf import Workflow
    from gwf.executors import Conda

    gwf = Workflow()

    gwf.target("Test", inputs=[], outputs=[], executor=Conda("myenv")) << """
    echo this will run inside the `myenv` Conda environment.
    """

  A default executor can also be specified for the workflow:

  .. code-block:: python

    from gwf import Workflow
    from gwf.executors import Conda

    gwf = Workflow(executor=Conda("myenv"))

    gwf.target("Test", inputs=[], outputs=[]) <<< """
    echo this will run inside the `myenv` Conda environment.
    """

  The available executors are documented :ref:`here <executors>`. By default,
  the ``Bash`` executor is used, so nothing will change for existing workflows.

  Executors are currently only available using the Slurm backend, but will be
  implemented for other backends in the near future.
* Writes to standard out and error inside the target run are now buffered. This
  prevents overloading network filesystems when a program is doing bad IO, such
  as printing a progress bar to standard error, which then has to be flushed to
  the filesystem very often. The buffer is written at least every 10 seconds,
  so progress indicators can still be used to check progress of a target in the
  log file, but it will not overload the filesystem.

  This functionality is also enabled by the new executor component.
* A backend for IBM Spectrum Load Sharing Facilility (LSF) was contributed by
  @gregoryleeman in PR #412 and #418.
* A backend for the Portable Batch System (PBS) was contributed in PR #419, also
  by  @gregoryleeman!
* You can now specify a `group` for a target which enables grouping targets:

  .. code-block:: python

    from gwf import Workflow

    gwf = Workflow()

    gwf.target("UnzipPatient1", ..., group="patient1") << "..."
    gwf.target("AnalyzePatient1", ..., group="patient1") << "..."
    gwf.target("UnzipPatient2", ..., group="patient2") << "..."
    gwf.target("AnalyzePatient2", ..., group="patient2") << "..."

  You can then tell ``gwf status`` to group the targets:

  .. code-block:: shell

    $ gwf status -f grouped
    patient1      0 shouldrun   1 submitted   0 running     0 completed   0 failed      1 cancelled
    patient2      0 shouldrun   0 submitted   0 running     0 completed   0 failed      2 cancelled

  This makes it easier to monitor large workflows. The feature was contributed
  by @jakobjn in PR #414 and extended by @dansondergaard.
* The ``gwf run`` and ``gwf status`` commands will now show an identifier that
  identifies each target in the backend. Concretely, this means that Slurm job
  id's are now easily obtainable:

  .. code-block:: shell

    $ gwf run SayHello
    Submitted target SayHello (id: 53839420)
    $ gwf status SayHello
    - SayHello      submitted (id: 53839420)

* Two new flags have been added to ``gwf run`` for controlling what is submitted.

  If ``--force`` is given, *gwf* will submit the specified targets and their
  dependencies no matter what.

  If ``--no-deps`` is given, *gwf* will only submit the specified targets, but
  not their dependencies. This can be combined with ``force`` to force submit
  a list of targets, but not any of their dependencies.
* Support for Python 3.13 has been added.

Changed
-------

* Support for Python 3.7 has been dropped.

Version 2.0.5
=============

Added
-----

* Added support for Python 3.12.

Version 2.0.4
=============

Fixed
-----

* Fixed bug introduced in 2.0.3 that only affected the Slurm backend.
* Fixed the behavior of the cancel command (can now cancel multiple target, even if cancelling one of them fails).

Version 2.0.3
=============

Fixed
-----

* Fixed path validation when using named paths (thanks, Anders!)

Changed
-------

* Updated list of known Slurm states in Slurm backend.

Version 2.0.2
=============

Added
-----

* You can now filter by "failed" status in the ``status`` command (thanks Ludvig).
* The output of ``gwf status -f summary`` is now more consistent, has colors,
  and shows all possible states instead of only non-zero ones (thanks, Jakob Grove).

Fixed
-----

* Fixed bug where a failed target would not be included as dependency when
  submitting an upstream target.
* When accounting is enabled in the Slurm backend (on by default), **gwf** will
  now requests job states from both ``squeue`` and ``sacct``, with the status
  from squeue taking precedence. This prevents strange-looking output from ``gwf
  status`` when a target has been submitted, but its state has not yet been
  flushed to the database by Slurm.

Version 2.0.1
=============

Added
-----

* You can now specify which parts of a workflow to touch with the `gwf touch` command.
  E.g. ``gwf touch Foo`` will only touch Foo and its dependencies (and their dependencies).

Fixed
-----

* Fixed loading of entrypoints on Python 3.7.

Version 2.0.0
=============

This is a major release. Many internals have been changed, specifically on the scheduler,
which is now much more robust and consistent, and to backends, which now can report whether
a target failed or not. The local backend has also been rewritten from scratch. In the
future it will be able to support e.g. target cancellation.

Added
-----

* For supported backends (currently `slurm` and `local`), the status command can now show
  failed targets. These will be rescheduled as if they were incomplete. If accounting is
  not enabled in slurm, this functionality must be turned off with
  `gwf config set backend.slurm.accounting_enabled no`.

Fixed
------

* The configuration directory will now always be placed next to the workflow file in use.
* The status command is now consistent with the output of `gwf run`n.
* The scheduler will not re-run targets that are not "should run".
* Targets with only outputs will now be scheduled if an output file is missing. Otherwise,
  it will not be scheduled.
* And many, many other tiny things!

Changed
-------

* Workflows can no longer be included into each other. All namespacing has been removed.
* Templates returning tuples are no longer supported. Templates must now return an ``AnonymousTarget``.
* APIs for loading workflow, using the scheduler, creating and loading backends etc. have
  been changed. This only affects plugins and custom backends.
* The `gwf init` command has been removed. Instead, *gwf* will ask if a skeleton should be
  created in the current directory, if it can't find an existing workflow file.
* Spec hashing is now turned off by default since it can be confusing. It can be turned on
  with `gwf config set use_spec_hashes yes`.

Version 1.8.5
=============

Added
-----

* Spec hashes will now be computed automatically on the first run with spec
  hashes enabled.
* ``gwf touch`` will now also update spec hashes.

Version 1.8.4
=============

Various tiny fixes to the scheduler, nothing user visible.

Version 1.8.3
=============

Fixed
-----

* Status output regression reported in #399 has been fixed.

Version 1.8.2
=============

Fixed
-----

* You can now protect an output file from being deleted by the ``gwf clean``
  command with the ``protect`` argument when creating a target. Specify a list of
  the files to be protected and ``gwf clean`` will ignore them (even ``gwf clean --all``).
  This was implemented before, but it was broken, so this is considered a fix.
* You can now load modules next to your workflow file, e.g., ``templates.py`` again.

Changed
-------

* Including one workflow into another workflow had very few real use cases, so the
  feature has been dropped.

Version 1.8.1
=============

Added
-----

* Will now recursively look for a workflow file in parent directories. This
  means that you can now run e.g. ``gwf status`` in any subdirectory below the
  directory containing your workflow file.

Fixed
-----

* Fixed bug in the scheduling code which caused workflows that ran just fine to
  report "shouldrun" for some targets.

Changed
-------

* The scheduler now outputs a data structure with all information from the
  scheduling process and with separate types for each "reason".
* The ``gwf status`` command is now more informative about the state of each
  target.

Version 1.8.0
=============

Changed
-------

* Drop support for Python <3.7.
* Use a nicer theme for the website!
* Removed the automatic update checking as it was fragile and usually not
  useful.
* If no backend is configured, **gwf** will not try to guess which backend to
  use based on which commands are available on the system.

Added
-----

* Pretty output for the ``info`` command with ``--format pretty``. Still
  supports the JSON output format as the default.
* Added ``clean_logs`` setting which defaults to `yes` (the current behaviour).
  Setting this to `no` will not remove old logs when ``gwf run`` is executed.
* Re-run targets when their spec changes (#180). Can be disabled with
  ``gwf config set use_spec_hashes no``. Enabled by default.

Fixed
-----

* Remove use of ``time.clock()`` as it was deprecated.
* Fix import of ``collections.abc.Mapping``.
* Replaced use of ``click.get_terminal_size`` with ``shutil.get_terminal_size``.


Version 1.7.2
=============

Fixed
-----

* Trying to cancel a target that is not running or completed will now *not* fail.
* Documented the Slurm target options ``mail_user`` and ``mail_type``.

Added
-----

* Slurm backend now accepts a ``gres`` target option that maps directly to
  Slurm's ``--gres`` flag.


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
