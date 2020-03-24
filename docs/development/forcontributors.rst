================
For Contributors
================

We appreciate all contributions to *gwf*, not just contributions to the code!
Think something is missing from the documentation? Defined useful snippets for
your text editor? Add it and submit a pull request!


Set Up a Development Environment
================================

We recommend that you set up a separate environment for *gwf* development, for
example with `venv <https://docs.python.org/3/library/venv.html>`_ or
`Conda <https://conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_.

Make Your Changes
=================

1. Fork the repository, clone it and create a branch for your changes::

    git checkout -b my-change

2. Make the necessary changes and add unit tests if necessary.

3. Add a description of the changes to ``CHANGELOG.rst`` and add yourself to
   ``CONTRIBUTORS.rst`` (if you're not already there).

4. Test your changes and check for style violations::

    make init      # to install gwf for development
    gwf ...        # test your changes by running gwf
    make lint      # to check for style issues
    make test      # to run tests
    make coverage  # to check test coverage

5. If everything is alright, commit your changes::

    git add .
    git commit -m "Added some-feature"


Show Us Your Contribution!
==========================

1. Push your committed changes back to your fork on GitHub::

    git push origin HEAD

2. Follow `these <https://help.github.com/articles/creating-a-pull-request/>`_ steps to create a pull request.

3. Check for comments and suggestions on your pull request and keep an eye on the
   `CI output <https://github.com/gwforg/gwf/actions?query=workflow%3A%22Run+tests%22>`_.
