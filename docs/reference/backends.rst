.. _backends:

Backends
========

*gwf* supports multiple backends for running workflows. If you don't find a
backend that suits your needs here, it's easy to write your own backend
<writing_backends>`. The source code for the built-in backends is a great
source of inspiration.

By default, *gwf* comes with the `local`, `slurm`, and `sge` backends.

Local
-----

.. automodule::
    gwf.backends.local


.. _slurm_backend:

Slurm
-----

.. automodule::
    gwf.backends.slurm


Sun Grid Engine (SGE)
---------------------

.. automodule::
    gwf.backends.sge
