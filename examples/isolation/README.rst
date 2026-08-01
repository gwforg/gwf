Target isolation
================

Run the workflow with::

    gwf run

``isolation=True`` enables the default copy-input/move-output configuration,
including blocked network access, as used by the two intentionally failing
targets below. The successful targets use explicit configurations to
demonstrate both input transfer modes. Set ``network=True`` on an
``IsolationConfig`` when a target needs normal network access.

``CopyInput`` copies ``inputs/samples/sample.txt`` into its isolated directory,
while ``SymlinkInput`` creates a symlink to
``inputs/references/reference.txt``. Their specifications verify the staging
method before reading the input.

The first two targets each create one declared and one undeclared output. Only
``results/copied/output.txt`` and ``results/symlinked/output.txt`` are
moved from the temporary directories to those paths in the workflow directory.
The files under ``work/copied`` and ``work/symlinked`` are discarded with the
temporary directories.

``MissingOutput`` intentionally demonstrates the failure case. It declares
``results/missing-output.txt`` but creates only
``work/missing-output/produced-but-undeclared.txt``. After its command exits
successfully, *gwf* detects the missing declared output and marks the target as
failed. The file that was created is not moved to the workflow directory, and
the declared output remains absent.

Nested file hierarchies
-----------------------

``NestedHierarchy`` declares two named inputs in separate subdirectories:
``inputs/samples/sample.txt`` and ``inputs/references/reference.txt``. *gwf*
recreates that hierarchy in the isolated directory, and also prepares the
parent directories for the declared ``results/nested/combined.txt`` output.
The target can therefore use the same relative paths in its specification.

The target also creates ``work/nested/intermediate.txt``. That entire ``work``
tree is discarded because the intermediate file is not a declared output.

Missing inputs
--------------

``MissingInput`` declares ``inputs/samples/does-not-exist.txt``. Since that
input does not exist and no other target produces it, *gwf* marks the target as
skipped before starting ``gwf.exec``. Its specification is never executed, so
neither ``results/missing-input.txt`` nor ``missing-input-command-ran.txt`` is
created.
