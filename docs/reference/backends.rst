.. _backends:

Backends
========

*gwf* supports multiple backends for running workflows. If you don't find a
backend that suits your needs here, it's easy to
:ref:`write your own backend <writing_backends>`.

By default, *gwf* comes with the `local` and `slurm` backends.

Local
-----

.. autoclass::
    gwf.backends.local.LocalBackend


Slurm
-----

.. autoclass::
   gwf.backends.slurm.SlurmBackend
