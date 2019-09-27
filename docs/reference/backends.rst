.. _backends:

Backends
========

*gwf* supports multiple backends for running workflows. If you don't find a
backend that suits your needs here, it's easy to
:ref:`write your own backend <writing_backends>`.

By default, *gwf* comes with the `local`, `slurm`, and `sge` backends.

Local
-----

.. autoclass::
    gwf.backends.local.LocalBackend


.. _slurm_backend:

Slurm
-----

.. autoclass::
    gwf.backends.slurm.SlurmBackend


Sun Grid Engine (SGE)
---------------------

.. autoclass::
    gwf.backends.sge.SGEBackend
