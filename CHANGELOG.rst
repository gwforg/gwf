Change Log
==========

Version 1.0 (Unreleased)
------------------------

First stable release of *gwf*! We strongly encourage users of pre-1.0 users to read the tutorial, since quite a lot
of things have changed. We also recommend reading the guide for converting pre-1.0 workflows to version 1.0. However,
users attempting to do this should be aware that the the template mechanism in 1.0 is slightly different and thus
requires rewriting template functions.


Version 1.0b10
--------------

Fixed
^^^^^

* Fixed a subtle bug in scheduling which caused problems when resubmitting a workflow where some targets were already running (a5d884).

* Fixed a bug in the ``SlurmBackend`` which caused *gwf* to crash if the Slurm queue contained a job with many dependencies (eb4446).

* Added back the `-e` flag in the ``logs`` command.

Version 1.0b9
-------------

Fixed
^^^^^

* Fixed a bug in the ``SlurmBackend`` which caused running targets as unknown (33a6bd).

Changed
^^^^^^^

* The Slurm backend's database of tracked jobs is now cleaned on initialization to keep it from growing indefinitely (bd3f95).

Version 1.0b8
-------------

Fixed
^^^^^

* Fixed a bug which caused the *gwf logs* command to always show stderr (01b267).

* Fixed a bug which caused dependencies to be set incorrectly when two targets depended on the same target (4d9e07).

Changed
^^^^^^^

* Improved error message when trying to create a target from an invalid template (d27d1f).

* Improved error message when assigning a non-string spec to a target (2aca0a).

* `gwf logs` command now outputs logs via a pager when the system supports it, unless `--no-pager` is used (01b267).

Added
^^^^^

* Added more tests to cover scenarios with included workflows when building the workflow graph (86a68d0).

* Added a bunch of documentation (69e136, 51a0e7, 942b05).

Version 1.0b7
-------------

Fixed
^^^^^

* Fixed bug in scheduling which was actually the cause of the incorrect scheduling that was "fixed" in 1.0b6.
  Also added documentation for ``gwf.core.schedule`` (7c47cb).

Changed
^^^^^^^

* Updated documentation in a bunch of places, mostly styling.

Version 1.0b6
-------------

Fixed
^^^^^

* A bug in ``SlurmBackend`` which caused dependencies between targets to not be set correctly (6b71d2).

Changed
^^^^^^^

* More improvements to and clean up of build process.
* Updated some examples in the tutorial with current output from *gwf* (42c5da).
* Logging output is now more consistent (b95af04).

Added
^^^^^

* Documentation for maintainers on how to merge in contributions and rolling a new release (fe1ee3).

Version 1.0b5
-------------

Fixed
^^^^^

* Unset option passed to backend causes error (#166, dcff44).
* Set import path to allow import of module in workflow file (64841c).

Changed
^^^^^^^

* Vastly improved build and deploy process. We're now actually building and testing with conda.
