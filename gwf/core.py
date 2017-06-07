import copy
import collections
import inspect
import itertools
import logging
import os
import os.path
import stat
import subprocess
import sys
from collections import defaultdict
from glob import glob as _glob
from glob import iglob as _iglob

from .backends.base import Status
from .exceptions import (CircularDependencyError,
                         FileProvidedByMultipleTargetsError,
                         FileRequiredButNotProvidedError, IncludeWorkflowError,
                         InvalidNameError, TargetExistsError, InvalidTypeError, InvalidPathError)
from .utils import (cache, load_workflow, is_valid_name, iter_inputs, iter_outputs, timer, parse_path,
                    match_targets)

logger = logging.getLogger(__name__)


def _norm_path(working_dir, path):
    path = str(path)
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(working_dir, path))


def _norm_paths(working_dir, paths):
    return [_norm_path(working_dir, path) for path in paths]


def _is_valid_list(obj):
    return isinstance(obj, collections.Sequence) and not isinstance(obj, str)


def normalized_paths_property(name):
    """Define a normalized paths property.

    This function can be used to define a property on an object which expects
    the value to be a list of paths. When the attribute is used, this function
    ensures that the paths will always be normalized and absolute.

    For an example of its usage, see the `inputs` and `outputs` attributes
    on :class:`~gwf.Target`. For a more general description of this approach
    to reusable attributes, see:

        http://chimera.labs.oreilly.com/books/1230000000393/ch09.html#_problem_164
    """
    storage_name = '_' + name

    @property
    def prop(self):
        return _norm_paths(self.working_dir, getattr(self, storage_name))

    @prop.setter
    def prop(self, value):
        setattr(self, storage_name, value)
    return prop


class Target(object):
    """Represents a target.

    A target is a unit of work that can be tied to other targets through
    dependencies on files.

    :ivar str name:
        Name of the target.
    :ivar list inputs:
        A list of input paths for this target.
    :ivar list outputs:
        A list of output paths for this target.
    :ivar dict options:
        Options such as number of cores, memory requirements etc. Options are
        backend-dependent. Backends will ignore unsupported options.
    :ivar str working_dir:
        Working directory of this target.
    :ivar str spec:
        The specification of the target.
    """

    inputs = normalized_paths_property('inputs')
    outputs = normalized_paths_property('outputs')

    def __init__(self, name, inputs, outputs, options, working_dir, spec=''):
        self.name = name
        if not is_valid_name(self.name):
            raise InvalidNameError('Target defined with invalid name: "{}".'.format(self.name))

        self.options = options
        self.working_dir = working_dir

        if not _is_valid_list(inputs):
            raise InvalidTypeError(
                'The argument `inputs` to target `{}` must be a list or tuple, not a string.'.format(name))
        if not _is_valid_list(outputs):
            raise InvalidTypeError(
                'The argument `outputs` to target `{}` must be a list or tuple, not a string.'.format(name))

        self.inputs = inputs
        self.outputs = outputs

        self.spec = spec

    def qualname(self, namespace):
        if namespace is not None:
            return '{}.{}'.format(namespace, self.name)
        return self.name

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
        if not isinstance(spec, str):
            msg = (
                'Target spec must be a string, not {}. '
                'Did you attempt to assign a template to this target? '
                'This is no is not allowed since version 1.0. Use the '
                'Workflow.target_from_template() method instead. See '
                'the tutorial for more details.'
            )
            raise InvalidTypeError(msg.format(type(spec)))
        self.spec = spec
        return self

    def __repr__(self):
        format_str = (
            '{}(name={!r}, inputs={!r}, outputs={!r}, options={!r}, '
            'working_dir={!r}, spec={!r})'
        )

        return format_str.format(
            self.__class__.__name__,
            self.name,
            self.inputs,
            self.outputs,
            self.options,
            self.working_dir,
            self.spec,
        )

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
            raise InvalidNameError('Workflow defined with invalid name: "{}".'.format(self.name))

        self.targets = {}
        self.defaults = defaults or {}

        self.working_dir = working_dir
        if self.working_dir is None:
            # Get the frame object of whatever called the Workflow.__init__
            # and extract the path of the file which is was defined in. Then
            # normalize the path and get the directory of the file.
            filename = inspect.getfile(sys._getframe(1))
            self.working_dir = os.path.dirname(os.path.realpath(filename))

    def _add_target(self, target, namespace=None):
        if namespace is not None:
            target.name = target.qualname(namespace)
        if target.name in self.targets:
            raise TargetExistsError(target)
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

        :param str name: Name of the target.
        :param tuple template: Target specification of the form (inputs, outputs, options, spec).

        Any further keyword arguments are passed to the backend and will
        override any options provided by the template.
        """
        try:
            inputs, outputs, template_options, spec = template
            if not isinstance(template_options, dict) or not isinstance(spec, str):
                raise Exception()
        except Exception:
            raise InvalidTypeError('Target `{}` received an invalid template.'.format(name))

        new_target = Target(
            name=name,
            inputs=inputs,
            outputs=outputs,
            options=options,
            working_dir=self.working_dir,
            spec=spec,
        )
        new_target.inherit_options(template_options)
        new_target.inherit_options(self.defaults)

        self._add_target(new_target)
        return new_target

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
            raise IncludeWorkflowError(
                'The included workflow has not been assigned a name. To '
                'include the workflow you must assign a name to the included '
                'workflow or set the namespace argument.'
            )
        namespace_prefix = namespace or other_workflow.name
        if namespace_prefix == self.name:
            raise IncludeWorkflowError('The included workflow has the same name as this workflow.')

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
            self.include_workflow(getattr(other_workflow, 'gwf'), namespace=namespace)
        else:
            raise TypeError('First argument must be either a string or a Workflow object.')

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
        return subprocess.check_output(*args, shell=True, cwd=self.working_dir, **kwargs)

    def __repr__(self):
        return '{}(name={!r}, working_dir={!r})'.format(
            self.__class__.__name__, self.name, self.working_dir
        )


class Graph(object):
    """Represents a finalized workflow graph.

    The represents the targets present in a workflow, but also their dependencies and the files they provide. When a
    graph is initialized it computes all dependency relations between targets, ensuring that the graph is semantically
    sane. Therefore, construction of the graph is an expensive operation which may raise a number of exceptions:

    :raises gwf.exceptions.FileProvidedByMultipleTargetsError:
        Raised if the same file is provided by multiple targets.
    :raises gwf.exceptions.FileRequiredButNotProvidedError:
        Raised if a target has an input file that does not exist on the
        file system and that is not provided by another target.
    :raises gwf.exceptions.InvalidPathError:
        Raised if a target has declared a directory in either `inputs` or `outputs`.
    :raises gwf.exceptions.CircularDependencyError:
        Raised if the workflow contains a circular dependency.

    If the graph is constructed successfully, the following instance variables will be available:

    :ivar dict targets:
        A dictionary mapping target names to instances of :class:`gwf.Target`.
    :ivar dict provides:
        A dictionary mapping a file path to the target that provides that path.
    :ivar dict dependencies:
        A dictionary mapping a target to a set of its dependencies.
    :ivar dict dependents:
        A dictionary mapping a target to a list of all targets which depend on the target.

    The graph can be manipulated in arbitrary, diabolic ways after it has been constructed. Checks are only performed
    at construction-time, thus introducing e.g. a circular dependency by manipulating *dependencies* will not raise an
    exception.
    """

    def __init__(self, targets):
        self.targets = targets

        logger.debug('Preparing workflow with %s targets defined.', len(self.targets))

        # The order is important here!
        self.provides = self.prepare_file_providers()
        self.dependencies = self.prepare_dependencies()
        self.dependents = self.prepare_dependents()

        self._check_for_circular_dependencies()

        self.file_cache = self.prepare_file_cache()
        logger.debug('Cached %d files.', len(self.file_cache))

    @timer('Prepared file providers in %.3fms.', logger=logger)
    def prepare_file_providers(self):
        provides = {}
        for target, path in iter_outputs(self.targets.values()):
            if path in provides:
                raise FileProvidedByMultipleTargetsError(
                    path, provides[path].name, target
                )

            provides[path] = target
        return provides

    @timer('Prepared dependencies in %.3fms.', logger=logger)
    def prepare_dependencies(self):
        dependencies = defaultdict(set)
        for target, path in iter_inputs(self.targets.values()):
            if path not in self.provides:
                if not os.path.exists(path):
                    raise FileRequiredButNotProvidedError(path, target)
                continue
            dependencies[target].add(self.provides[path])
        return dependencies

    @timer('Prepared dependents in %.3fms.', logger=logger)
    def prepare_dependents(self):
        dependents = defaultdict(list)
        for target, deps in self.dependencies.items():
            for dep in deps:
                dependents[dep].append(target)
        return dependents

    @timer('Prepared file cache in %.3fms.', logger=logger)
    def prepare_file_cache(self):
        cache = {}

        input_iter = iter_inputs(self.targets.values())
        output_iter = iter_outputs(self.targets.values())
        for target, path in itertools.chain(input_iter, output_iter):
            if path in cache:
                continue

            try:
                st = os.stat(path)
            except FileNotFoundError:
                cache[path] = None
            else:
                if stat.S_ISDIR(st.st_mode):
                    raise InvalidPathError('Path {} used in {} points to a directory.'.format(target.name, path))
                cache[path] = st.st_mtime
        return cache

    @timer('Checked for circular dependencies in %.3fms.', logger=logger)
    def _check_for_circular_dependencies(self):
        """Check for circular dependencies in the graph.

        Raises :class:`CircularDependencyError` if a circular dependency is found.
        """
        fresh, started, done = 0, 1, 2

        nodes = self.targets.values()
        state = dict((n, fresh) for n in nodes)

        def visitor(node):
            state[node] = started
            for dep in self.dependencies[node]:
                if state[dep] == started:
                    raise CircularDependencyError(node)
                elif state[dep] == fresh:
                    visitor(dep)
            state[node] = done

        for node in nodes:
            if state[node] == fresh:
                visitor(node)

    @cache
    def should_run(self, target):
        """Return whether a target should be run or not."""
        if any(self.should_run(dep) for dep in self.dependencies[target]):
            logger.debug(
                '%s should run because one of its dependencies should run.',
                target.name
            )
            return True

        if target.is_sink:
            logger.debug(
                '%s should run because it is a sink.',
                target.name
            )
            return True

        if any(self.file_cache[path] is None for path in target.outputs):
            logger.debug(
                '%s should run because one of its output files does not exist.',
                target.name
            )
            return True

        if target.is_source:
            logger.debug(
                '%s should not run because it is a source.',
                target.name
            )
            return False

        youngest_in_ts, youngest_in_path = max(
            (self.file_cache[path], path) for path in target.inputs
        )

        logger.debug(
            '%s is the youngest input file of %s with timestamp %s.',
            youngest_in_path,
            target.name,
            youngest_in_ts
        )

        oldest_out_ts, oldest_out_path = min(
            (self.file_cache[path], path) for path in target.outputs
        )

        logger.debug(
            '%s is the oldest output file of %s with timestamp %s.',
            oldest_out_path,
            target.name,
            youngest_in_ts
        )

        if youngest_in_ts > oldest_out_ts:
            logger.debug(
                '%s should run since %s is larger than %s.',
                target.name,
                youngest_in_ts,
                oldest_out_ts
            )
            return True

        return False

    def endpoints(self):
        """Return a set of all targets that are not depended on by other targets."""
        return set(self.targets.values()) - set(self.dependents.keys())

    def find(self, names):
        """Return a set of targets matching the given names/patterns."""
        return match_targets(names, self.targets)

    def __iter__(self):
        return iter(self.targets.values())


def schedule(graph, backend, target, dry_run=False):
    """Schedule a target and its dependencies.

    Scheduling a target will determine whether the target needs to run based on
    whether it already has been submitted and whether any of its dependencies have
    been submitted.

    Targets that should run will be submitted to *backend*, unless *dry_run*
    is set to ``True``.

    Returns ``True`` if *target* was submitted to the backend (even when
    *dry_run* is ``True``).

    :param gwf.Graph graph:
        Graph of the workflow.
    :param gwf.backends.Backend backend:
        An instance of :class:`gwf.backends.Backend` to which targets will be
        submitted.
    :param gwf.Target target:
        Target to be scheduled.
    :param bool dry_run:
        If ``True``, targets will not be submitted to the backend. Defaults
        to ``False``.
    """
    logger.info('Scheduling target %s.', target.name)

    if backend.status(target) != Status.UNKNOWN:
        logger.debug('Target %s has already been submitted.', target.name)
        return True

    submitted_deps = set()
    for dependency in sorted(graph.dependencies[target], key=lambda t: t.name):
        was_submitted = schedule(graph, backend, dependency, dry_run=dry_run)
        if was_submitted:
            submitted_deps.add(dependency)

    if submitted_deps or graph.should_run(target):
        if dry_run:
            logger.info('Would submit target %s.', target.name)
        else:
            logger.info('Submitting target %s.', target.name)
            backend.submit(target, dependencies=submitted_deps)
        return True
    else:
        logger.debug('Target %s should not run.', target.name)
        return False


def schedule_many(graph, backend, targets, **kwargs):
    """Schedule multiple targets and their dependencies.

    This is a convenience function for scheduling multiple targets. See
    :func:`schedule` for a detailed description of the arguments and behavior.

    :param list targets: A list of targets to be scheduled.
    """
    schedules = []
    for target in targets:
        was_submitted = schedule(graph, backend, target, **kwargs)
        schedules.append(was_submitted)
    return schedules
