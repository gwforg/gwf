.. _configuration:

=============
Configuration
=============

Configuration of *gwf* is project-specific and thus all configuration must be done in the project directory where the
workflow file is located.

To see the value of a configuration key, use:

.. code-block:: bash

    $ gwf config get KEY

To set the value of a key (or update it, if it already exists):

.. code-block:: bash

    $ gwf config set KEY VALUE

Note that a keys are often of the form `this.is.a.key`. For example, the local backend supports the `local.port` setting
which sets the port that the workers are running on. To set this settings, just run:

.. code-block:: bash

    $ gwf config set local.port 4321

Now, when you run *gwf* with the local backend, it will try to connect workers on port 4321.

Your configuration is stored in the current working directory, which will usually be your project directory, in a file
called ``.gwfconf.json``. This means that all configuration is project-specific, which helps with reproducibility. You
can inspect and change the file directly, but this is not recommended unless you really know what you're doing.

Core settings are listed on the :ref:`settings` page. To see which options are available for a specific backend, refer
to the :ref:`backends` documentation.
