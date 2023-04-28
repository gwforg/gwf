==============
Welcome to gwf
==============

*gwf* is a flexible, pragmatic workflow tool for building and running large,
scientific workflows. It runs on Python 3.7+ and is developed at GenomeDK,
Aarhus University.

Examples
  To get a feeling for what a *gwf* workflow looks like, have a look at a few
  `examples <https://github.com/gwforg/gwf/tree/master/examples>`_.

Getting started
  To quickly get started writing workflows in *gwf* you can read the
  :ref:`tutorial`.

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

Installation
============

We'll assume that you have the Anaconda_ distribution installed and that you are
familiar with how to install and manage packages and environments through the
*conda* package manager.

To install *gwf* via Conda::

    conda config --add channels gwforg
    conda install gwf

We recommend that you install *gwf* in a project-specific environment::

    conda config --add channels gwforg
    conda create -n myproject gwf dep1 dep2 ...
    conda activate myproject

You can find the code for *gwf* `here <https://github.com/gwforg/gwf>`_. You are
encouraged to report any issues through the
`issue tracker <https://github.com/gwforg/gwf/issues>`_, which is also a good place
to ask questions.

.. _Anaconda: https://www.continuum.io/downloads

User's Guide
============

.. toctree::
  :maxdepth: 3

  guide/index

Extending
=========

.. toctree::
  :maxdepth: 3

  plugins
  reference/index

Release history
===============

.. toctree::
  :maxdepth: 2

  changelog
