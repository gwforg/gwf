Welcome to gwf
==============

|status| |docs| |license| |versions| |format| |ci-status| |coveralls| |downloads|

*gwf* is a flexible, pragmatic workflow tool for building and running large,
scientific workflows developed at the Bioinformatics Research Centre
(BiRC), Aarhus University.

*gwf* provides a make-like system for specifying dependencies between tasks to
be run. Dependencies are based on files, but unlike *make* where several input
files can be used for dependencies for one output file, the workflow supports
tasks that takes several input files and produces several output files. Tasks
then have dependencies on selected files, either produced by one or many other
tasks or assumed to be present on the system.

To quickly get started with writing workflows in *gwf* you can read the
:ref:`tutorial`. If you want to extend *gwf* you should read
:ref:`writing_plugins`. We don't have the backend you need to run your workflow
on your cluster? See the :ref:`writing_backends` section!

Installation
------------

To install *gwf* via conda::

    conda install -c dansondergaard gwf

We recommend that you install *gwf* in a project-specific environment::

    conda create -n myproject python=3.5 dep1 dep2 ...
    source activate myproject
    conda install -c dansondergaard gwf

You can find the code for *gwf* `here <https://github.com/mailund/gwf>`_. You are
encouraged to report any issues through the
`issue tracker <https://github.com/mailund/gwf/issues>`_.

Upgrading from an Older Version
-------------------------------

If you have been using an older version of *gwf* you will have to
make a few changes to your workflows to be able to run them with version 1.0
and above. :doc:`Read more... </upgrading_from_pre_1.0>`

Getting Started
---------------

.. toctree::
   :maxdepth: 2

   tutorial
   configuration
   best_practices
   backends
   plugins

Extending gwf
-------------

.. toctree::
   :maxdepth: 2

   writingplugins
   writingbackends
   contributing

API Reference
-------------

.. toctree::
   :maxdepth: 2

   api


.. |ci-status| image:: 	https://img.shields.io/travis/mailund/gwf.svg
    :target: https://travis-ci.org/mailund/gwf
    :alt: Build status
.. |docs| image:: https://readthedocs.org/projects/gwf/badge/?version=latest&style=flat
    :target: http://gwf.readthedocs.io
    :alt: Documentation
.. |format| image:: https://img.shields.io/pypi/format/gwf.svg
    :target: https://pypi.python.org/pypi/gwf
    :alt: Kit format
.. |downloads| image:: https://img.shields.io/pypi/dm/gwf.svg
    :target: https://pypi.python.org/pypi/gwf
    :alt: Monthly PyPI downloads
.. |versions| image:: https://img.shields.io/pypi/pyversions/gwf.svg
    :target: https://pypi.python.org/pypi/gwf
    :alt: Python versions supported
.. |status| image:: https://img.shields.io/pypi/status/gwf.svg
    :target: https://pypi.python.org/pypi/gwf
    :alt: Package stability
.. |license| image:: https://img.shields.io/pypi/mailund/gwf.svg
    :target: https://pypi.python.org/pypi/gwf
    :alt: License
.. |coveralls| image:: https://img.shields.io/coveralls/mailund/gwf.svg
    :target: https://coveralls.io/github/mailund/gwf
    :alt: Coverage
