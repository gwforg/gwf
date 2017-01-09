import inspect
import itertools
import logging
import os.path
import sys
from collections import defaultdict

from .exceptions import (CircularDependencyError,
                         FileProvidedByMultipleTargetsError,
                         FileRequiredButNotProvidedError, IncludeWorkflowError,
                         InvalidNameError, TargetExistsError)
from .utils import (cache, dfs, get_file_timestamp, import_object,
                    is_valid_name, iter_inputs, iter_outputs, merge, timer)

logger = logging.getLogger(__name__)


def _norm_path(working_dir, path):
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(working_dir, path))


def _norm_paths(working_dir, paths):
    return [_norm_path(working_dir, path) for path in paths]


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
        backend-dependent.
    :ivar gwf.Workflow workflow:
        Workflow that the target is associated to.
    :ivar str spec:
        The specification of the target.
    """

    inputs = normalized_paths_property('inputs')
    outputs = normalized_paths_property('outputs')

    def __init__(self, name, inputs, outputs, options, workflow, namespace=None, spec=''):
        self.name = name
        if not is_valid_name(self.name):
            raise InvalidNameError(
                'Target defined with invalid name: "{}".'.format(self.name)
            )

        self.options = options
        self.workflow = workflow

        self.inputs = inputs
        self.outputs = outputs

        self.spec = spec

    def qualname(self, namespace):
        if namespace is not None:
            return '{}.{}'.format(namespace, self.name)
        if self.workflow.name is not None:
            return '{}.{}'.format(self.workflow.name, self.name)
        return self.name

    @property
    def working_dir(self):
        return self.workflow.working_dir

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
            self.inputs = options.pop('inputs', list)
            self.outputs = options.pop('outputs', list)
            self.options = options
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

    def target(self, name, inputs=None, outputs=None, **options):
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

        if inputs is None:
            inputs = []
        if outputs is None:
            outputs = []

        new_target = Target(
            name, inputs, outputs, options, workflow=self
        )

        self._add_target(new_target)
        return new_target

    def include_path(self, path, namespace=None):
        """Include targets from another :class:`gwf.Workflow` into this workflow.

        See :func:`~gwf.Workflow.include`.
        """
        other_workflow = import_object(path)
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
            self._add_target(target, namespace)

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

    def __repr__(self):
        return '{}(name={!r}, working_dir={!r}, targets={!r})'.format(
            self.__class__.__name__, self.name, self.working_dir, self.targets
        )


class PreparedWorkflow(object):

    """Represents a finalized workflow graph.

    If :class:`gwf.PreparedWorkflow` is initialized with the *workflow*
    parameter, the :class:`gwf.PreparedWorkflow` calls :meth:`prepare` with the
    workflow.

    :ivar targets: initial value: dict()
    :ivar working_dir: initial value: None
    :ivar provides: initial value: None
    :ivar dependencies: initial value: None
    :ivar dependents: initial value: None
    """

    def __init__(self, workflow=None, config=None, backend=None):
        self.targets = {}
        self.working_dir = None

        self.provides = None
        self.dependencies = None
        self.dependents = None
        self.file_cache = None

        if workflow is not None:
            self.prepare(workflow, config, backend)

    def prepare(self, workflow, config, backend):
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
        self.workflow = workflow
        self.config = config or {}
        self.backend = backend

        self.targets = workflow.targets
        self.working_dir = workflow.working_dir

        logger.debug(
            'Preparing workflow with %s targets defined.',
            len(self.targets)
        )

        logger.debug('Received configuration: %r.', config)

        # The order is important here!
        self.provides = self.prepare_file_providers()
        self.dependencies = self.prepare_dependencies()
        self.dependents = self.prepare_dependents()

        self._check_for_circular_dependencies()
        self._inherit_target_options()

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
        for target in self.targets.values():
            for dep in self.dependencies[target]:
                if target in dfs(dep, self.dependencies):
                    raise CircularDependencyError(target)

    def _inherit_target_options(self):
        for target_name, target in self.targets.items():
            self.targets[target_name].options = merge(self.config, {
                option: value
                for option, value in target.options.items()
            })

            self.targets[target_name].options = {
                option: value
                for option, value in target.options.items()
                if option in self.backend.supported_options
            }

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

    def __repr__(self):
        return '{}(working_dir={!r}, targets={!r})'.format(
            self.__class__.__name__, self.working_dir, self.targets
        )


class Event:
    """A simple event implementation.

    This class implements the observer pattern in a reusable fashion. Event
    objects can be created anywhere (in modules, class definitions etc.) and
    anyone can register for the event or trigger it.

    .. note::
        Even though it is possible to trigger any event, please only trigger
        your own events unless you really know what you're doing.

    Triggering an event will call all registered callbacks with the arguments
    given to :func:`trigger`.
    """

    def __init__(self, name):
        self.name = name
        self._callbacks = set()
        self._logger = logging.getLogger('{}.{}'.format(__name__, self.name))

    def register(self, callback):
        """Register a callback to this event.

        :param callable callback:
            a callable do call whenever the event is triggered.

        Registering the same callback twice will only add the callback once.
        """
        if callback not in self._callbacks:
            self._logger.debug('Registered callback: %r.', callback)
            self._callbacks.add(callback)

    def unregister(self, callback):
        """Unregister a callback.

        :param callable callback:
            a callable to unregister.

        It is safe to unregister the same callback twice. The second call to
        :func:`unregister` will be a no-op.
        """
        self._logger.debug('Unregistered callback: %r.', callback)
        self._callbacks.discard(callback)

    def trigger(self, *args, **kwargs):
        """Trigger the event.

        All arguments will be forwarded to the callbacks. Callbacks are called
        in no particular order.
        """
        for callback in self._callbacks:
            self._logger.debug(
                'Triggering callback: %r with args %r and kwargs %r.',
                callback,
                args,
                kwargs,
            )
            callback(*args, **kwargs)

    def __repr__(self):
        return '{}(name={})'.format(self.__class__.__name__, self.name)
