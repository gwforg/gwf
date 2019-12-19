import copy
import collections
import collections.abc
import inspect
import subprocess
import sys
import warnings
import os.path
import unicodedata
from glob import glob as _glob
from glob import iglob as _iglob

from .compat import fspath
from .exceptions import WorkflowError, NameError, TypeError
from .utils import is_valid_name, parse_path, load_workflow


def _has_nonprintable_char(s):
    chars = enumerate((unicodedata.category(char) == "Cc", char) for char in s)
    for pos, (unprintable, char) in chars:
        if unprintable:
            return (
                s.encode("unicode_escape").decode("utf-8"),
                char.encode("unicode_escape").decode("utf-8"),
                pos,
            )
    return None


def _check_path(path, target_name, mode):
    if not path:
        msg = 'Target "{}" has an empty {} path.'.format(target_name, mode)
        raise WorkflowError(msg)

    result = _has_nonprintable_char(path)
    if result is not None:
        clean_path, char, pos = result
        msg = (
            'Path "{}" in target "{}" {}s contains a '
            'non-printable character "{}" on position {}. '
            "This is always unintentional and can cause "
            "strange behaviour."
        ).format(clean_path, target_name, mode, char, pos)
        raise WorkflowError(msg)
    return path


def _norm_path(working_dir, path):
    path = fspath(path)
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(working_dir, path))


def _norm_paths(working_dir, paths):
    return [_norm_path(working_dir, path) for path in paths]


def _flatten(t):
    res = []

    def flatten_rec(g):
        if isinstance(g, str) or hasattr(g, "__fspath__"):
            res.append(g)
        elif isinstance(g, collections.abc.Mapping):
            for k, v in g.items():
                flatten_rec(v)
        else:
            for v in g:
                flatten_rec(v)

    flatten_rec(t)
    return res


def select(lst, fields):
    """Select fields from an iterable of dictionaries.

    Given an iterable of dictionaries and an iterable of key names, selects the
    given key names from the dictionaries and returns an iterable of
    dictionaries containing only those keys.

    For example, given the list::

        >>> lst = [
        ...     {'A': 'a1', 'B': 'b1'},
        ...     {'A': 'a2', 'B': 'b2'},
        ...     {'A': 'a3', 'B': 'b3'},
        ... ]

    We can select only the `A` keys from the dictionaries with::

        >>> list(select(lst, ['A']))
        [{'A': 'a1'}, {'A': 'a2'}, {'A': 'a3'}]

    :param iterable lst:
        An iterable of dictionaries.
    :param iterable fields:
        An iterable of key names.
    """
    fields = tuple(fields)
    for item in lst:
        dct = {}
        for name in fields:
            dct[name] = item[name]
        yield dct


def collect(lst, fields, rename=None):
    """Collect values from an iterable of dictionaries into a dictionary.

    Given an iterable of dictionaries and an iterable of key names, collect the
    value of each field name in `fields` and return a dictionary of lists for
    each field.

    For example, given the list:

        >>> lst = [
        ...     {'A': 'a1', 'B': 'b1'},
        ...     {'A': 'a2', 'B': 'b2'},
        ...     {'A': 'a3', 'B': 'b3'},
        ... ]

    We can collect the values into a dictionary of lists with::

        >>> collect(lst, ['A'])
        {'As': ['a1', 'a2', 'a3']}

    :param iterable lst:
        An iterable of dictionaries.
    :param iterable fields:
        An iterable of key names.
    """
    fields = tuple(fields)
    if rename is None:
        rename = {}
    selected = collections.defaultdict(list)
    for item in lst:
        for name in fields:
            selected_name = rename.get(name, name + "s")
            selected[selected_name].append(item[name])
    return dict(selected)


class TargetList(list):
    """A list of target objects with access to all inputs and outputs.

    This is a thin wrapper around a normal list and thus provides all normal
    ``list`` methods. However, it provides access to the collective
    inputs and outputs of the targets contained in the list.
    """

    @property
    def outputs(self):
        """Return a list of the outputs of all targets.

        The returned list may be a list of strings, lists or dictionaries
        depending on the form of the outputs of the contained targets.
        """
        return [target.outputs for target in self]

    @property
    def inputs(self):
        """Return a list of the inputs of all targets.

        The returned list may be a list of strings, lists or dictionaries
        depending on the form of the inputs of the contained targets.
        """
        return [target.inputs for target in self]

    def __str__(self):
        class_name = self.__class__.__name__
        if not self:
            return "{}(targets=[])".format(class_name)
        return "{}(targets=[{!r}, ...])".format(class_name, self[0])

    def __repr__(self):
        return str(self)


class AnonymousTarget:
    """Represents an unnamed target.

    An anonymous target is an unnamed, abstract target much like the tuple
    returned by function templates. Thus, `AnonymousTarget` can also be used as
    the return value of a template function.

    :ivar list inputs:
        A string, list or dictionary containing inputs to the target.
    :ivar list outputs:
        A string, list or dictionary containing outputs to the target.
    :ivar dict options:
        Options such as number of cores, memory requirements etc. Options are
        backend-dependent. Backends will ignore unsupported options.
    :ivar str working_dir:
        Working directory of this target.
    :ivar str spec:
        The specification of the target.
    :ivar set protect:
        An iterable of protected files which will not be removed during
        cleaning, even if this target is not an endpoint.
    """

    _creation_order = 0

    def __init__(
        self, inputs, outputs, options, working_dir=None, spec="", protect=None
    ):
        self.options = options
        self.working_dir = working_dir
        self.inputs = inputs
        self.outputs = outputs
        self._spec = spec

        self.order = AnonymousTarget._creation_order
        AnonymousTarget._creation_order += 1

        if protect is None:
            self.protected = set()
        else:
            self.protected = set(protect)

    @property
    def spec(self):
        return self._spec

    @spec.setter
    def spec(self, value):
        if not isinstance(value, str):
            msg = (
                "Target spec must be a string, not {}. Did you attempt to "
                "assign a template to this target? This is no is not allowed "
                "since version 1.0. Use the Workflow.target_from_template() "
                "method instead. See the tutorial for more details."
            )
            raise TypeError(msg.format(type(value)))

        self._spec = value

    @property
    def is_source(self):
        """Return whether this target is a source.

        A target is a source if it does not depend on any files.
        """
        return not self.inputs

    @property
    def is_sink(self):
        """Return whether this target is a sink.

        A target is a sink if it does not output any files.
        """
        return not self.outputs

    def inherit_options(self, super_options):
        options = super_options.copy()
        options.update(self.options)
        self.options = options

    def __lshift__(self, spec):
        self.spec = spec
        return self

    def __repr__(self):
        return "{}(inputs={!r}, outputs={!r}, options={!r}, working_dir={!r}, spec={!r})".format(
            self.__class__.__name__,
            self.inputs,
            self.outputs,
            self.options,
            self.working_dir,
            self.spec,
        )

    def __str__(self):
        return "{}_{}".format(self.__class__.__name__, id(self))


class Target(AnonymousTarget):
    """Represents a target.

    This class inherits from :class:`AnonymousTarget`.

    A target is a named unit of work that declare their file *inputs* and
    *outputs*. Target names must be valid Python identifiers.

    A script (or spec) is associated with the target. The script must be a
    valid Bash script and should produce the files declared as *outputs* and
    consume the files declared as *inputs*. Both parameters must be provided
    explicitly, even if no inputs or outputs are needed. In that case, provide
    the empty list::

        Target('Foo', inputs=[], outputs=[], options={}, working_dir='/tmp')

    The *inputs* and *outputs* arguments can either be a string, a list or
    a dictionary. If a dictionary is given, the keys act as names for the
    files. The values may be either strings or a list of strings::

        foo = Target(
            name='foo',
            inputs={'A': ['a1', 'a2'], 'B': 'b'},
            outputs={'C': ['a1b', 'a2b], 'D': 'd},
        )

    This is useful for referring the outputs of a target::

        bar = Target(
            name='bar',
            inputs=foo.outputs['C'],
            outputs='result',
        )

    The target can also specify an *options* dictionary specifying the
    resources needed to run the target. The options are consumed by the backend
    and may be ignored if the backend doesn't support a given option. For
    example, we can set the *cores* option to set the number of cores that the
    target uses::

        Target('Foo', inputs=[], outputs=[], options={'cores': 16}, working_dir='/tmp')

    To see which options are supported by your backend of choice, see the
    documentation for the backend.

    :ivar str name:
        Name of the target.

    .. versionchanged:: 1.6.0
        Named inputs and outputs were added. Prior versions require *inputs*
        and *outputs* to be lists.
    """

    def __init__(
        self, name, inputs, outputs, options, working_dir=None, spec="", protect=None
    ):
        self.name = name
        if not is_valid_name(self.name):
            raise NameError('Target defined with invalid name: "{}".'.format(self.name))

        _check_path(working_dir, target_name=self.name, mode="working_dir")
        for path in inputs:
            _check_path(path, target_name=self.name, mode="input")
        for path in outputs:
            _check_path(path, target_name=self.name, mode="output")

        super().__init__(
            inputs=inputs,
            outputs=outputs,
            options=options,
            working_dir=working_dir,
            spec=spec,
            protect=protect,
        )

    def flattened_inputs(self):
        return _norm_paths(self.working_dir, _flatten(self.inputs))

    def flattened_outputs(self):
        return _norm_paths(self.working_dir, _flatten(self.outputs))

    @classmethod
    def empty(cls, name):
        """Return a target with no inputs, outputs and options.

        This is mostly useful for testing.
        """
        return cls(
            name=name, inputs=[], outputs=[], options={}, working_dir=os.getcwd()
        )

    def qualname(self, namespace):
        if namespace is not None:
            return "{}.{}".format(namespace, self.name)
        return self.name

    def __repr__(self):
        return "{}(name={!r}, ...)".format(self.__class__.__name__, self.name)

    def __str__(self):
        return self.name


class Workflow(object):
    """Represents a workflow.

    This is the most central user-facing abstraction in *gwf*.

    A workflow consists of a collection of targets and has methods for adding
    targets to the workflow in two different ways. A workflow can be initialized
    with the following arguments:

    :ivar str name: initial value: None
        The name is used for namespacing when including workflows. See
        :func:`~include` for more details on namespacing.
    :ivar str working_dir:
        The directory containing the file where the workflow was initialized.
        All file paths used in targets added to this workflow are relative to
        the working directory.
    :ivar dict defaults:
        A dictionary with defaults for target options.

    By default, *working_dir* is set to the directory of the workflow file which
    initialized the workflow. However, advanced users may wish to set it manually.
    Targets added to the workflow will inherit the workflow working directory.

    The *defaults* argument is a dictionary of option defaults for targets and
    overrides defaults provided by the backend. Targets can override the
    defaults individually. For example::

        gwf = Workflow(defaults={
            'cores': 12,
            'memory': '16g',
        })

        gwf.target('Foo', inputs=[], outputs=[]) << \"\"\"echo hello\"\"\"
        gwf.target('Bar', inputs=[], outputs=[], cores=2) << \"\"\"echo world\"\"\"

    In this case `Foo` and `Bar` inherit the `cores` and `memory` options set in
    `defaults`, but `Bar` overrides the `cores` option.

    See :func:`~include` for a description of the use of the `name` argument.
    """

    def __init__(self, name=None, working_dir=None, defaults=None):
        self.name = name
        if self.name is not None and not is_valid_name(self.name):
            raise NameError(
                'Workflow defined with invalid name: "{}".'.format(self.name)
            )

        self.targets = {}
        self.defaults = defaults or {}

        self.working_dir = working_dir
        if self.working_dir is None:
            # Get the frame object of whatever called the Workflow.__init__
            # and extract the path of the file which is was defined in. Then
            # normalize the path and get the directory of the file.
            filename = inspect.getfile(sys._getframe(1))
            self.working_dir = os.path.dirname(os.path.realpath(filename))

    @classmethod
    def from_path(cls, path):
        """Return workflow object for the workflow given by `path`.

        Returns a :class:`~gwf.Workflow` object containing the workflow object
        of the workflow given by `path`.

        :arg str path: Path to a workflow file, optionally specifying a
            workflow object in that file.
        """
        basedir, filename, obj = parse_path(path)
        return load_workflow(basedir, filename, obj)

    @classmethod
    def from_config(cls, config):
        """Return workflow object for the workflow specified by `config`.

        See :func:`Workflow.from_path` for further information.
        """
        return cls.from_path(config["file"])

    def _add_target(self, target, namespace=None):
        if namespace is not None:
            target.name = target.qualname(namespace)
        if target.name in self.targets:
            raise WorkflowError(
                'Target "{}" already exists in workflow.'.format(target.name)
            )
        self.targets[target.name] = target

    def target(self, name, inputs, outputs, **options):
        """Create a target and add it to the :class:`gwf.Workflow`.

        This is syntactic sugar for creating a new :class:`~gwf.Target` and
        adding it to the workflow. The target is also returned from the method
        so that the user can directly manipulate it, if necessary. For example,
        this allows assigning a spec to a target directly after defining it::

            workflow = Workflow()
            workflow.target('NewTarget', inputs=['test.txt', 'out.txt']) <<< '''
            cat test.txt > out.txt
            echo hello world >> out.txt
            '''

        This will create a new target named `NewTarget`, add it to the workflow
        and assign a spec to the target.

        :param str name: Name of the target.
        :param iterable inputs: List of files that this target depends on.
        :param iterable outputs: List of files that this target produces.

        Any further keyword arguments are passed to the backend.
        """
        new_target = Target(
            name=name,
            inputs=inputs,
            outputs=outputs,
            options=options,
            working_dir=self.working_dir,
        )
        new_target.inherit_options(self.defaults)

        self._add_target(new_target)
        return new_target

    def target_from_template(self, name, template, **options):
        """Create a target from a template and add it to the :class:`gwf.Workflow`.

        This is syntactic sugar for creating a new :class:`~gwf.Target` and
        adding it to the workflow. The target is also returned from the method
        so that the user can directly manipulate it, if necessary.

        .. code-block:: python

            workflow = Workflow()
            workflow.target_from_template('NewTarget', my_template())

        This will create a new target named `NewTarget`, configure it based
        on the specification in the template `my_template`, and
        add it to the workflow.

        :param str name:
            Name of the target.
        :param AnonymousTarget template:
            The anonymous target which describes the template.

        Any further keyword arguments are passed to the backend and will
        override any options provided by the template.
        """
        if isinstance(template, AnonymousTarget):
            new_target = Target(
                name=name,
                inputs=template.inputs,
                outputs=template.outputs,
                options=options,
                working_dir=template.working_dir or self.working_dir,
                spec=template.spec,
            )

            new_target.inherit_options(template.options)
        elif isinstance(template, tuple):
            warnings.warn(
                (
                    "Creating a target from a tuple template is deprecated, "
                    "and will be removed in gwf 2.0. Make your template "
                    "function return an AnonymousTarget instead."
                ),
                DeprecationWarning,
            )

            try:
                inputs, outputs, template_options, spec = template
            except ValueError:
                raise TypeError(
                    "Target `{}` received an invalid template.".format(name)
                )

            new_target = Target(
                name=name,
                inputs=inputs,
                outputs=outputs,
                options=options,
                working_dir=self.working_dir,
                spec=spec,
            )

            new_target.inherit_options(template_options)
        else:
            raise TypeError("Target `{}` received an invalid template.".format(name))

        new_target.inherit_options(self.defaults)
        self._add_target(new_target)
        return new_target

    def map(self, template_func, inputs, extra=None, name=None, **kwargs):
        """Add targets to the workflow given a template and a list of inputs.

        This method accepts a template function and an iterable of inputs. For
        each item in `inputs` it produces a target using the template function
        and adds the target to this workflow.

        For example, given this template:

        .. code-block::

            def copy_file(from_file):
                inputs = {'from_file': from_file}
                outputs = {'to_file': to_file + '.copy'}
                options = {}
                spec = "cp {from_file} {to_file}".format(from_file, to_file)
                return AnonymousTarget(
                    inputs=inputs,
                    outputs=outputs,
                    options=options,
                    spec=spec
                )

        and this list of files:

        .. code-block::

            files = ['file1', 'file2', 'file3']

        we can generate targets to copy all three files:

        .. code-block::

            gwf = Workflow()
            res = gwf.map(copy_file, files)

        The :func:`map` method returns a :class:`TargetList` which contains the
        generated targets.

        :param template_func:
            A function or callable class instance that returns an
            :class:`AnonymousTarget`. Essentially a *template function*.
        :param iterable inputs:
            An iterable of inputs for the generated targets. This can be an
            iterable of strings, tuples or dictionaries.
        :param mapping extra:
            A mapping of extra keyword arguments to be passed to the template.
        :param name: Must be either `None`, a string or a function.

            If `None` is given, the name of each target will be generated from
            the name of the template and an index.

            If a string is given, e.g. `foo`, the generated names will be
            `foo_0`, `foo_1`, etc.

            If a function is given, it must have the signature
            `f(idx, target)` where `idx` is the index and `target` is the
            :class:`AnonymousTarget` returned by the template. The function
            must return the name to assign to the target as a string.

        Any remaining keyword arguments will be passed directly to
        :func:`target_from_template` and thus override template-specified
        target options.
        """

        if extra is None:
            extra = {}

        if not (
            callable(template_func)
            and (
                hasattr(template_func, "__name__")
                or hasattr(template_func, "__class__")
            )
        ):
            raise ValueError(
                "Argument `template_func` must be a function or a callable class instance."
            )

        def template_namer(idx, target):
            if hasattr(template_func, "__name__"):
                name = template_func.__name__
            else:
                name = template_func.__class__.__name__
            return "{name}_{idx}".format(name=name, idx=idx)

        def string_namer(idx, target):
            return "{name}_{idx}".format(name=name, idx=idx)

        if name is None:
            name_func = template_namer
        elif isinstance(name, str):
            name_func = string_namer
        else:
            name_func = name

        targets = TargetList()
        for idx, args in enumerate(inputs):
            if isinstance(args, collections.abc.Mapping):
                template = template_func(**args, **extra)
            elif isinstance(args, collections.abc.Iterable) and not isinstance(
                args, str
            ):
                template = template_func(*args, **extra)
            else:
                template = template_func(args, **extra)

            target_name = name_func(idx, template)
            target = self.target_from_template(
                name=target_name, template=template, **kwargs
            )
            targets.append(target)
        return targets

    def include_path(self, path, namespace=None):
        """Include targets from another :class:`gwf.Workflow` into this workflow.

        See :func:`~gwf.Workflow.include`.
        """
        basedir, filename, obj = parse_path(path)
        other_workflow = load_workflow(basedir, filename, obj)
        self.include_workflow(other_workflow, namespace=namespace)

    def include_workflow(self, other_workflow, namespace=None):
        """Include targets from another :class:`gwf.Workflow` into this workflow.

        See :func:`~gwf.Workflow.include`.
        """
        if other_workflow.name is None and namespace is None:
            raise WorkflowError(
                "The included workflow has not been assigned a name. To "
                "include the workflow you must assign a name to the included "
                "workflow or set the namespace argument."
            )
        namespace_prefix = namespace or other_workflow.name
        if namespace_prefix == self.name:
            raise WorkflowError(
                "The included workflow has the same name as this workflow."
            )

        for target in other_workflow.targets.values():
            self._add_target(copy.deepcopy(target), namespace=namespace_prefix)

    def include(self, other_workflow, namespace=None):
        """Include targets from another :class:`gwf.Workflow` into this workflow.

        This method can be given either an :class:`gwf.Workflow` instance,
        a module or a path to a workflow file.

        If a module or path the workflow object to include will be determined
        according to the following rules:

        1. If a module object is given, the module must define an attribute
           named `gwf` containing a :class:`gwf.Workflow` object.
        2. If a path is given it must point to a file defining a module with an
           attribute named `gwf` containing a :class:`gwf.Workflow`
           object. If you want to include a workflow with another name you can
           specify the attribute name with a colon, e.g.::

                /some/path/workflow.py:myworkflow

           This will include all targets from the workflow `myworkflow`
           declared in the file `/some/path/workflow.py`.

        When a :class:`gwf.Workflow` instance has been obtained, all
        targets will be included directly into this workflow. To avoid name
        clashes the `namespace` argument must be provided. For example::

            workflow1 = Workflow()
            workflow1.target('TestTarget')

            workflow2 = Workflow()
            workflow2.target('TestTarget')

            workflow1.include(workflow2, namespace='wf1')

        The workflow now contains two targets named `TestTarget` (defined in
        `workflow2`) and `wf1.TestTarget` (defined in `workflow1`). The
        `namespace` parameter can be left out if the workflow to be included
        has been named::

            workflow1 = Workflow(name='wf1')
            workflow1.target('TestTarget')

            workflow2 = Workflow()
            workflow2.target('TestTarget')

            workflow1.include(workflow2)

        This yields the same result as before. The `namespace` argument can be
        used to override the specified name::

            workflow1 = Workflow(name='wf1')
            workflow1.target('TestTarget')

            workflow2 = Workflow()
            workflow2.target('TestTarget')

            workflow1.include(workflow2, namespace='foo')

        The workflow will now contain targets named `TestTarget` and
        `foo.TestTarget`.
        """
        if isinstance(other_workflow, Workflow):
            self.include_workflow(other_workflow, namespace=namespace)
        elif isinstance(other_workflow, str):
            self.include_path(other_workflow, namespace=namespace)
        elif inspect.ismodule(other_workflow):
            self.include_workflow(getattr(other_workflow, "gwf"), namespace=namespace)
        else:
            raise TypeError(
                "First argument must be either a string or a Workflow object."
            )

    def glob(self, pathname, *args, **kwargs):
        """Return a list of paths matching `pathname`.

        This method is equivalent to :func:`python:glob.glob`, but searches with
        relative paths will be performed relative to the working directory
        of the workflow.
        """
        if not os.path.isabs(pathname):
            pathname = os.path.join(self.working_dir, pathname)
        return _glob(pathname, *args, **kwargs)

    def iglob(self, pathname, *args, **kwargs):
        """Return an iterator which yields paths matching `pathname`.

        This method is equivalent to :func:`python:glob.iglob`, but searches with
        relative paths will be performed relative to the working directory
        of the workflow.
        """
        if not os.path.isabs(pathname):
            pathname = os.path.join(self.working_dir, pathname)
        return _iglob(pathname, *args, **kwargs)

    def shell(self, *args, **kwargs):
        """Return the output of a shell command.

        This method is equivalent to :func:`python:subprocess.check_output`, but
        automatically runs the command in a shell with the current working
        directory set to the working directory of the workflow.

        .. versionchanged:: 1.0

            This function no longer return a list of lines in the output, but a
            byte array with the output, exactly like :func:`python:subprocess.check_output`.
            You may specifically set *universal_newlines* to `True` to get a
            string with the output instead.
        """
        return subprocess.check_output(
            *args, shell=True, cwd=self.working_dir, **kwargs
        )

    def __repr__(self):
        return "{}(name={!r}, working_dir={!r})".format(
            self.__class__.__name__, self.name, self.working_dir
        )
