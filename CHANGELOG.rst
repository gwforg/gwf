Change Log
==========

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

* A bug in `SlurmBackend` which caused dependencies between targets to not be set correctly (6b71d2).

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
