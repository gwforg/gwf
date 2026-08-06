.. _target_isolation:

==========================
Isolating Target Execution
==========================

.. versionchanged:: 3.0.0
    Target isolation was added.

Target isolation gives a target a deliberately constrained execution
environment. Its purpose is to make workflows more reproducible by exposing
hidden dependencies: a target should depend on its declared inputs, create its
declared outputs, and only use network access when that access is intentional.

.. caution::
    Target isolation is *not* enabled by default.

For an isolated target, *gwf* constructs a working directory from the declared
inputs and blocks network access by default. After a successful run, only
declared outputs are transferred back to the workflow directory, while other
files created by the target are discarded.

Enable isolation for a target by passing an :class:`~gwf.IsolationConfig`:

.. code-block:: python

    from gwf import IsolationConfig

    gwf.target(
        "Analyze",
        inputs=["data/input.txt"],
        outputs=["results/output.txt"],
        isolation=True,
    ) << """
    analyze data/input.txt > results/output.txt
    """

For the default isolation settings, ``isolation=True`` is equivalent to
``isolation=IsolationConfig()``, which disables network access and runs the
target in an isolated, temporary directory with input files symlinked and
declared outputs moved back to the workflow directory after a successful run.

If ``isolation`` is ``None`` (the default), the target runs in its original
working directory without isolation.

The isolated directory exists only for the duration of the target. *gwf*
removes it after the command finishes and the declared outputs have been
handled. Set ``network=True`` to retain normal network access for a target.

Inputs and outputs must currently be regular files. Directory inputs and
outputs are not supported.

Declared outputs are copied or moved from the temporary directory to their
configured paths in the workflow directory only after the target command
completes successfully.

Input and output paths must be relative paths within the workflow directory.
Absolute paths and paths containing a ``..`` component are not supported for
an isolated target. Commands should use the same relative paths that are
declared in ``inputs`` and ``outputs``.

Configuring File Transfers
==========================

The transfer method can be configured for each target. The valid options are:

* ``inputs="symlink"`` creates a symlink to each declared input (default).
* ``inputs="copy"`` copies each declared input into the temporary directory (
  also known as staging).
* ``outputs="move"`` moves each declared output to its path in the workflow
  directory (default).
* ``outputs="copy"`` copies each declared output to its path in the workflow
  directory.

For example:

.. code-block:: python

    gwf.target(
        "Analyze",
        inputs=["data/input.txt"],
        outputs=["results/output.txt"],
        isolation=IsolationConfig(inputs="symlink", outputs="copy"),
    ) << """
    analyze data/input.txt > results/output.txt
    """

Copying inputs gives the target independent files. Symlinking is faster and
uses less temporary disk space, but the links still refer to the original
files.

Moving outputs avoids an extra copy and is therefore the default. Thus, the
default configuration is:

.. code-block:: python

    isolation=IsolationConfig(inputs="symlink", outputs="move", network=False)

While copying or symlinking inputs and copying or moving declared outputs,
*gwf* logs each operation and its duration to standard error. These messages
make it possible to identify slow file transfers and see which paths were
handled.

Choosing the Temporary Directory Root
=====================================

By default, isolated directories are created in the system temporary directory.
Set ``root`` to place them under another directory, such as node-local scratch
storage:

.. code-block:: python

    isolation=IsolationConfig(root="/scratch")

Each target gets a unique temporary directory below this root. Setting this to a
node-local scratch directory and setting ``inputs="copy"`` effectively enables
staging of the input files, which can help avoid a lot of I/O-related issues.

The configured root must already exist and be writable on the machine where the
target runs. *gwf* removes the target-specific directory after execution, but
does not remove the configured root.

A relative ``root`` is resolved relative to the target's original working
directory.

Workflow Defaults
=================

Isolation can be enabled for every target by setting a workflow default:

.. code-block:: python

    gwf = Workflow(defaults={
        "isolation": True,
    })

Individual targets may replace the default configuration:

.. code-block:: python

    gwf.target(
        "LargeInput",
        inputs=["large-input.dat"],
        outputs=["result.dat"],
        isolation=IsolationConfig(inputs="symlink", outputs="move"),
    ) << """
    process large-input.dat > result.dat
    """

Network Access
==============

.. caution::
    Blocking network access requires unprivileged user namespaces to be enabled
    on the execution host. This is commonly available on current Linux systems,
    but may be disabled by system administrators.

Network access is blocked by default for isolated targets. This restriction
also applies while an executor such as Conda, Pixi, Apptainer, or Singularity
starts the target. Environments and container images must therefore already be
available locally.

To allow normal network access for a target, set ``network=True``:

.. code-block:: python

    gwf.target(
        "Download",
        inputs=[],
        outputs=["downloaded.dat"],
        isolation=IsolationConfig(network=True),
    ) << """
    download-data > downloaded.dat
    """

Blocking network access by default provides some protection against supply-chain
attacks and malicious software, but be aware that this isolation is not a
full-blown security boundary.
