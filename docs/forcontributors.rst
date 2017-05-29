For contributors
================

We appreciate all contributions to *gwf*, not just contributions to the code! Think something is missing from the
documentation? Defined useful snippets for your text editor? Add it and submit a pull request!


Set up a development environment
--------------------------------

We strongly recommend that you use the `Anaconda Python distribution <https://www.continuum.io/anaconda-overview>`_
and the conda package manager to set up a development environment (actually, we recommend it for all of your Python
work).

1. Download and install the Anaconda Python distribution following the instructions
   `here <https://www.continuum.io/downloads>`_.

2. Create an environment for *gwf* development::

    conda create -n gwfdev python=3.4

3. Activate the environment::

    source activate gwfdev


Make your changes
-----------------

1. Fork the repository, clone it and create a branch for your changes::

    git checkout -b my-change

2. Make the necessary changes and add unit tests if necessary.

3. Add a description of the changes to ``CHANGELOG.rst`` and add yourself to ``CONTRIBUTORS.rst`` (if you're not
   already there).

4. Test your changes and check for style violations::

    make init      # to install gwf for development
    gwf ...        # test your changes by running gwf
    make lint      # to check for style issues
    make test      # to run tests
    make coverage  # to check test coverage

5. If everything is alright, commit your changes::

    git add .
    git commit -m "Added some-feature"


Show us your contribution!
--------------------------

1. Push your committed changes back to your fork on GitHub::

    git push origin HEAD

2. Follow `these <https://help.github.com/articles/creating-a-pull-request/>`_ steps to create a pull request.

3. Check for comments and suggestions on your pull request and keep an eye on the
   `CI output <https://travis-ci.org/mailund/gwf>`_.
