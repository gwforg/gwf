Upgrading from Pre-1.0
======================


Pre-1.0 versions of *gwf* used a slightly different syntax for defining
workflows. However, the changes are quite minimal and thus updating an
existing workflow to run with version 1.0 and up should be straightforward.
First, you should change the import from::

    from gwf import *

to::

    from gwf import Workflow, template

Right after this, insert this line::

    gwf = Workflow()

Now search for all occurrences of `target(` and replace them with `gwf.target(`.
Your workflow should now run as before. For example, a workflow that previously
looked like this::

    from gwf import *

    target('Foo', input=['hello.txt'], output=['bye.txt']) << """
    ...
    """

should now look like this::

    from gwf import Workflow, template

    gwf = Workflow()

    gwf.target('Foo', input=['hello.txt'], output=['bye.txt']) << """
    ...
    """
