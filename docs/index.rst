==============
Welcome to gwf
==============

*gwf* is a flexible, pragmatic workflow tool for building and running large,
scientific workflows. It runs on Python 3.5+ and is developed at the Bioinformatics
Research Centre (BiRC), Aarhus University.

Examples
  To get a feeling for what a *gwf* workflow looks like, have a look at a few
  `examples <https://github.com/mailund/gwf/tree/master/examples>`_.

Getting started
  To quickly get started writing workflows in *gwf* you can read the
  :ref:`tutorial`.

Extending
  We don't have the backend you need to run your workflow on your cluster?
  See the :ref:`writing_backends` section to roll your own.

Contributing
  We aim to make *gwf* a community developed project. Learn how to
  :doc:`contribute </development/forcontributors>`.

Installation
============

To install *gwf* via conda::

    conda config --add channels gwforg
    conda install gwf

We recommend that you install *gwf* in a project-specific environment::

    conda config --add channels gwforg
    conda create -n myproject python=3.5 gwf dep1 dep2 ...
    source activate myproject

You can find the code for *gwf* `here <https://github.com/gwforg/gwf>`_. You are
encouraged to report any issues through the
`issue tracker <https://github.com/gwforg/gwf/issues>`_, which is also a good place
to ask questions.

Getting Started
===============

.. toctree::
   :maxdepth: 2

   tutorial
   patterns

Topic Guides
============

.. toctree::
   :maxdepth: 2

   topic-guides/index

Development
===========

.. toctree::
   :maxdepth: 2

   development/index

Reference
=========

.. toctree::
   :maxdepth: 2

   reference/settings
   reference/backends
   reference/api

.. include:: ../CHANGELOG.rst
.. include:: ../CONTRIBUTORS.rst
