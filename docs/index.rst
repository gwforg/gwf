Welcome to gwf
==============

*gwf* is a flexible, pragmatic workflow tool for building and running large,
scientific workflows. It is developed at the Bioinformatics Research Centre
(BiRC), Aarhus University.

In *gwf*, a workflow consists of a collection of *targets*. A target declares
its input and output files

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
  :doc:`contribute </contributing>`.

Installation
------------

To install *gwf* via conda::

    conda install -c gwforg gwf

We recommend that you install *gwf* in a project-specific environment::

    conda create -n myproject python=3.5 dep1 dep2 ...
    source activate myproject
    conda install -c gwforg gwf

You can find the code for *gwf* `here <https://github.com/gwforg/gwf>`_. You are
encouraged to report any issues through the
`issue tracker <https://github.com/gwforg/gwf/issues>`_, which is also a good place
to ask questions.

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
   backends
   plugins

Extending gwf
-------------

.. toctree::
   :maxdepth: 2

   writingbackends
   contributing

API
---

.. toctree::
   :maxdepth: 2

   api
