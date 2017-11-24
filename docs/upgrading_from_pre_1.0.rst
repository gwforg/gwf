======================
Upgrading from Pre-1.0
======================


.. note::
    *gwf* only runs on Python 3.5+. This means that you may have to upgrade your
    Python installation. It is still possible to run workflows using scripts
    written in Python 2.7 using e.g. environments.


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
For example, a workflow that previously looked like this::

    from gwf import *

    target('Foo', input='hello.txt', output=['bye.txt']) << """
    ...
    """

should now look like this::

    from gwf import Workflow, template

    gwf = Workflow()

    gwf.target('Foo', input='hello.txt', output=['bye.txt']) << """
    ...
    """

In target definitions you should rename `input` to `inputs` and
`output` to `outputs`. The workflow then ends up looking like this::

    from gwf import Workflow, template

    gwf = Workflow()

    gwf.target('Foo', inputs='hello.txt', outputs=['bye.txt']) << """
    ...
    """

Finally, in older versions of *gwf* you were allowed to just pass a string as input
or output (as for `inputs` above). In 1.0, you must always pass an iterable. Thus,
the workflow ends up like this::

    from gwf import Workflow, template

    gwf = Workflow()

    gwf.target('Foo', inputs=['hello.txt'], outputs=['bye.txt']) << """
    ...
    """

Old versions of *gwf* provided a ``template()`` function for creating templates
and doing string substitution, as well as the ability to have functions return 
templates. The ``template()`` function has been removed, and template functions
are now the only way to create templates.

Additionally, function templates must now return the inputs and outputs of the
template separate from the other options. To see how, take a look at
:ref:`templates`.