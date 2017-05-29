Change Log
==========

Version 1.0b7
-------------

* FIXED: Fixed bug in scheduling which was actually the cause of the incorrect scheduling that was "fixed" in 1.0b6.
  Also added documentation for ``gwf.core.schedule`` (7c47cb).
* CHANGED: Updated documentation in a bunch of places, mostly styling.

Version 1.0b6
-------------

* FIXED: A bug in `SlurmBackend` which caused dependencies between targets to not be set correctly (6b71d2).
* CHANGED: More improvements to and clean up of build process.
* CHANGED: Updated some examples in the tutorial with current output from *gwf* (42c5da).
* CHANGED: Logging output is now more consistent (b95af04).
* ADDED: Documentation for maintainers on how to merge in contributions and rolling a new release (fe1ee3).

Version 1.0b5
-------------

* FIXED: #166 (dcff44).
* FIXED: import module in workflow file (64841c).
* CHANGED: Vastly improved build and deploy process. We're now actually building and testing with conda.
