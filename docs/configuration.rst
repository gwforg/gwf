.. _configuration:

Configuration
=============

Configuration of *gwf* is project-specific and thus all configuration must be done
in the project directory where the workflow file is located.

To see the value of a configuration key, use:

.. code-block:: bash

    $ gwf config get KEY

To set the value of a key (or update it, if it already exists):

.. code-block:: bash

    $ gwf config set KEY VALUE

Note that a keys are often of the form `this.is.a.key`.

Some backends can be configured in various ways. To see which options are available for
a specific backend, refer to the :ref:`backends` documentation.

Reference
---------

* **backend (str):** Set the backend. Corresponds to the ``--backend`` flag (default: local).
* **verbose (str):** Set the verbosity. Corresponds to the ``--verbose`` flag (default: info).
* **no_color (bool):** If `true`, colors will not be used (default: false).
