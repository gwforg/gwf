==============
Welcome to gwf
==============

*gwf* is a flexible, pragmatic workflow tool for building and running large,
scientific workflows. It runs on Python 3.5+ and is developed at the Bioinformatics
Research Centre (BiRC), Aarhus University.

Examples
  To get a feeling for what a *gwf* workflow looks like, have a look at a few
  `examples <https://github.com/gwforg/gwf/tree/master/examples>`_.

Getting started
  To quickly get started writing workflows in *gwf* you can read the
  :ref:`tutorial`.

Extending
  We don't have the backend you need to run your workflow on your cluster?
  See the :ref:`writing_backends` section to roll your own.

Contributing
  We aim to make *gwf* a community developed project. Learn how to
  :doc:`contribute </development/forcontributors>`.

Features
========

* Easy to adopt, there's no special syntax to learn if you already know Python
* Automatically resolves dependencies between targets based on filenames
* Only submits targets when their output files are not up to date
* Supports multiple backends like Slurm, SGE and a local backend for testing
* Fire-and-forget, does not require you to use *screen* or *tmux* to keep your
  workflow running
* Commands for cleaning temporary data from your workflow
* Friendly to your system administrator!

User's Guide
============

.. toctree::
  :maxdepth: 3

  guide/index

Development
===========

.. toctree::
  :maxdepth: 2

  development/index
  changelog

API Reference
=============

.. toctree::
  :maxdepth: 3

  reference/index


.. include:: ../CONTRIBUTORS.rst
