.. _formaintainers:

===============
For Maintainers
===============

The *gwf* build, testing and deployment process is automated through Travis.

Merging Changes
===============

1. Make sure that the changes have proper test coverage, e.g. by checking the branch
   on `Coveralls <https://coveralls.io/github/gwforg/gwf>`_.

2. Check that the PR includes necessary updates of :file:`CHANGELOG.rst` and
   :file:`CONTRIBUTORS.rst`.

3. Always make a merge commit (don't rebase/fast-forward). The merge commit will be
   referenced in the change log.

4. Add the change to the change log for the coming (draft) release on
   `GitHub <https://github.com/gwforg/gwf/releases>`_. Make sure to follow the
   formatting used in previous change logs. Also, read about
   `how to keep a change log <http://keepachangelog.com/en/0.3.0/>`_.

Rolling a New Release
=====================

1. Make sure that all changes for the new release have been merged into ``master``
   and that tests pass. Check the `CI output <https://github.com/gwforg/gwf/actions?query=workflow%3A%22Run+tests%22>`_.

2. Make any other release-related changes such as adding new contributors to
   ``CONTRIBUTORS.rst`` or adding missing items to ``CHANGELOG.rst``.

3. Increase the version number in ``gwf/__init__.py``

4. Commit the changes and push the branch. Wait for tests to run.

5. Make a new release by tagging the merge commit with the version number, e.g.
   ``vX.X.X``. Push the tag and wait for Travis to catch up.

6. Run ``make package``, then ``make publish`` to publish the source
   distribution and wheel to PyPI.

The documentation will be automatically be built.
