import inspect
import itertools
import logging
import os.path
import subprocess
import sys
from collections import defaultdict
from glob import glob as _glob
from glob import iglob as _iglob

import collections

from .exceptions import (CircularDependencyError,
                         FileProvidedByMultipleTargetsError,
                         FileRequiredButNotProvidedError, IncludeWorkflowError,
                         InvalidNameError, TargetExistsError, InvalidTypeError, TargetDoesNotExistError)
from .utils import (cache, dfs, get_file_timestamp, load_workflow,
                    is_valid_name, iter_inputs, iter_outputs, timer, parse_path)

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
    :ivar dict options:
        Options such as number of cores, memory requirements etc. Options are
        backend-dependent.
    :ivar gwf.Workflow workflow:
        Workflow that the target is associated to.
    :ivar str spec:
        The specification of the target.
    """

    inputs = normalized_paths_property('inputs')
    outputs = normalized_paths_property('outputs')

    def __init__(self, name, options, working_dir, spec=''):
        self.name = name
        if not is_valid_name(self.name):
            raise InvalidNameError('Target defined with invalid name: "{}".'.format(self.name))

        # we need these for __repr__ to work -- it sets the propertis
        # the actual inputs and outputs won't be set until we validate
        # the target -- we keep them in the options to make sure
        # that they are explicitly set.
        self.inputs = []
        self.outputs = []

        self.options = options.copy()
        self.working_dir = working_dir

        self.spec = spec

    def validate(self):
        """Validate that the target is correctly specified."""

        if 'inputs' not in self.options:
            raise InvalidTypeError(
                'Target `{}` does not specify its `inputs`.'.format(self.name))
        if 'outputs' not in self.options:
            raise InvalidTypeError(
                'Target `{}` does not specify its `outputs`.'.format(self.name))

        inputs = self.options['inputs']
        outputs = self.options['outputs']

        # after we get the inputs and outputs options we need to remove
        # them from the table so they don't get passed to backends.
        # FIXME: This is not ideal, since it means that you cannot validate
        # a target twice or you will get an error :(
        del self.options['inputs']
        del self.options['outputs']

        if not _is_valid_list(inputs):
            raise InvalidTypeError(
                'The argument `inputs` to target `{}` must be a list or tuple, not a string.'.format(self.name))
        if not _is_valid_list(outputs):
            raise InvalidTypeError(
                'The argument `outputs` to target `{}` must be a list or tuple, not a string.'.format(self.name))

        self.inputs = inputs
        self.outputs = outputs

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

    def __lshift__(self, spec):
        if isinstance(spec, tuple):
            options, spec = spec
            self.inputs = options.get('inputs', list)
            self.outputs = options.get('outputs', list)

            # Override template options with target options.
            self.working_dir = options.get('working_dir', self.working_dir)

            target_options = self.options.copy()
            self.options = options.copy()
            self.options.update(target_options)

            self.spec = spec
        else:
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

    :ivar str name: initial value: None
        The name is used for namespacing when including workflows. See
        :func:`~include` for more details on namespacing.
    :ivar dict targets:
        A dictionary of the targets in this workflow.
    :ivar str working_dir:
        The directory containing the file where the workflow was initialized.
        All file paths used in targets added to this workflow are relative to
        the working directory.
    """

    def __init__(self, name=None, working_dir=None, defaults=None):
        self.name = name
        if self.name is not None and not is_valid_name(self.name):
            raise InvalidNameError(
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

    def _add_target(self, target, namespace=None):
        qualname = target.qualname(namespace)
        if qualname in self.targets:
            raise TargetExistsError(target)
        self.targets[qualname] = target

    def target(self, name, **options):
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

        Any further keyword arguments are passed to the backend.
        """
        merged_options = self.defaults.copy()
        merged_options.update(options)

        new_target = Target(
            name, merged_options,
            working_dir=self.working_dir,
        )

        self._add_target(new_target)
        return new_target

    def include_path(self, path, namespace=None):
        """Include targets from another :class:`gwf.Workflow` into this workflow.

        See :func:`~gwf.Workflow.include`.
        """
        basedir, filename, obj = parse_path(path)
        other_workflow = load_workflow(basedir, filename, obj)
        self.include_workflow(other_workflow, namespace=namespace)
        return other_workflow

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
            raise IncludeWorkflowError(
                'The included workflow has the same name as this workflow.'
            )

        for target in other_workflow.targets.values():
            self._add_target(target, namespace_prefix)

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
            self.include_workflow(
                getattr(other_workflow, 'gwf'), namespace=namespace)
        else:
            raise TypeError('First argument must be either a string or a '
                            'Workflow object.')

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

        .. note:: This function has changed in 1.0. It will no longer return a
          list of lines in the output, but a byte array with the output,
          exactly like :func:`python:subprocess.check_output`. You may specifically
          set *universal_newlines* to `True` to get a string with the output
          instead.
        """
        return subprocess.check_output(*args, shell=True, cwd=self.working_dir)

    def __repr__(self):
        return '{}(name={!r}, working_dir={!r})'.format(
            self.__class__.__name__, self.name, self.working_dir
        )


class Graph(object):

    """Represents a finalized workflow graph.

    If :class:`gwf.Graph` is initialized with the *workflow*
    parameter, the :class:`gwf.Graph` calls :meth:`prepare` with the
    workflow.

    :ivar targets: initial value: dict()
    :ivar provides: initial value: None
    :ivar dependencies: initial value: None
    :ivar dependents: initial value: None
    """

    def __init__(self, targets):
        self.targets = targets
        for target in self.targets.values():
            target.validate()

        self.provides = None
        self.dependencies = None
        self.dependents = None
        self.file_cache = None

        self.prepare()

    def prepare(self):
        """Prepare this workflow given a :class:`gwf.Workflow` instance.

        :param gwf.Workflow workflow:
            The workflow which should be prepared.
        :raises FileProvidedByMultipleTargetsError:
            Raised if the same file is provided by multiple targets.
        :raises FileRequiredButNotProvidedError:
            Raised if a target has an input file that does not exist on the
            file system and that is not provided by another target.
        :raises CircularDependencyError:
            Raised if the workflow contains a circular dependency.
        """

        logger.debug(
            'Preparing workflow with %s targets defined.',
            len(self.targets)
        )

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
        dependencies = defaultdict(list)
        for target, path in iter_inputs(self.targets.values()):
            if path not in self.provides:
                if not os.path.exists(path):  # pragma: no branch
                    raise FileRequiredButNotProvidedError(path, target)
                continue  # pragma: no cover
            dependencies[target].append(self.provides[path])
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
        input_iter = iter_inputs(self.targets.values())
        output_iter = iter_outputs(self.targets.values())
        return {path: get_file_timestamp(path)
                for _, path in itertools.chain(input_iter, output_iter)}

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

    def dfs(self, target):
        """Return a depth-first search path starting at `target`."""
        return dfs(target, self.dependencies)

    def endpoints(self):
        """Return a set of all targets that are not depended on by other targets."""
        return set(self.targets.values()) - set(self.dependents.keys())

    def get_targets_by_name(self, names):
        matched_targets = []
        for name in names:
            if name not in self.targets:
                raise TargetDoesNotExistError(name)
            matched_targets.append(self.targets[name])
        return matched_targets


def schedule(graph, backend, target, dry_run=False):
    """Schedule and submit a :class:`gwf.Target` and its dependencies.

    This method is provided by :class:`Backend` and should not be overriden.
    """
    logger.info('Scheduling target %s.', target.name)

    if backend.submitted(target):
        logger.debug('Target %s has already been submitted.', target.name)
        return []

    scheduled = []
    for dependency in graph.dfs(target):
        if dependency.name != target.name:
            logger.info(
                'Scheduling dependency %s of %s.',
                dependency.name,
                target.name
            )

        if backend.submitted(dependency):
            logger.debug(
                'Target %s has already been submitted.',
                dependency.name
            )
            continue

        if not graph.should_run(dependency):
            logger.debug(
                'Target %s should not run.',
                dependency.name
            )
            continue

        if dry_run:
            logger.info('Would submit target %s.', dependency.name)
        else:
            logger.info('Submitting target %s.', dependency.name)
            backend.submit(dependency, graph.dependencies[dependency])
        scheduled.append(dependency)

    return scheduled


def schedule_many(graph, backend, targets, **kwargs):
    """Schedule a list of :class:`gwf.Target` and their dependencies.

    Will schedule the targets in `targets` with :func:`schedule`
    and return a list of schedules.

    This method is provided by :class:`Backend` and should not be overriden.

    :param list targets: A list of targets to be scheduled.
    :return: A list of schedules, one for each target in `targets`.
    """
    schedules = []
    for target in targets:
        schedules.append(schedule(graph, backend, target, **kwargs))
    return schedules
