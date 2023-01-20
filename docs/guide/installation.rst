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
    source activate myproject

You can find the code for *gwf* `here <https://github.com/gwforg/gwf>`_. You are
encouraged to report any issues through the
`issue tracker <https://github.com/gwforg/gwf/issues>`_, which is also a good place
to ask questions.

.. _Anaconda: https://www.continuum.io/downloads
