.. _configuration:

=============
Configuration
=============

Configuration of *gwf* is project-specific and thus all configuration must be
done in the project directory where the workflow file is located.

Commands for Configuration
==========================

To see the value of a configuration key, use:

.. code-block:: bash

    $ gwf config get KEY

To set the value of a key (or update it, if it already exists):

.. code-block:: bash

    $ gwf config set KEY VALUE

Note that a keys are often of the form `this.is.a.key`. For example, the local
backend supports the `local.port` setting which sets the port that the workers
are running on. To set this settings, just run:

.. code-block:: bash

    $ gwf config set local.port 4321

Now, when you run *gwf* with the local backend, it will try to connect workers
on port 4321.

Your configuration is stored in the current working directory, which will
usually be your project directory, in a file called ``.gwfconf.json``. This
means that all configuration is project-specific, which helps with
reproducibility. You can inspect and change the file directly, but this is not
recommended unless you really know what you're doing.

Core settings are listed in the section below. To see which options are
available for a specific backend, refer to the :ref:`backends` documentation.

Available Settings
==================

This page lists settings that are used by *gwf*. Backends and plugins may define
their own settings, but these are documented for each backend/plugin
individually. See :ref:`configuration` if in doubt about how to configure *gwf*.

* **backend (str):** Set the backend. Corresponds to the ``--backend`` flag (default: `local`).
* **verbose (str):** Set the verbosity. Corresponds to the ``--verbose`` flag (default: `info`).
* **no_color (bool):** If `true`, colors will not be used (default: `false`).
* **check_updates (bool):** When `true`, *gwf* will automatically check for
  updates at most once per day. Setting this to `false` will deactivate the
  check completely.
