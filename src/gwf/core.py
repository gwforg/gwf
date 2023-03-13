import hashlib
import json
import logging
import os
import os.path
import unicodedata
from collections import defaultdict
from collections.abc import Mapping
from enum import Enum
from os import fspath

import attrs

from .exceptions import GWFError, NameError, WorkflowError
from .utils import is_valid_name, timer

logger = logging.getLogger(__name__)


def _flatten(t):
    res = []

    def flatten_rec(g):
        if isinstance(g, str) or hasattr(g, "__fspath__"):
            res.append(g)
        elif isinstance(g, Mapping):
            for k, v in g.items():
                flatten_rec(v)
        else:
            for v in g:
                flatten_rec(v)

    flatten_rec(t)
    return res


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


def hash_spec(spec):
    return hashlib.sha1(spec.encode("utf-8")).hexdigest()


def hash_specs(graph):
    hashes = {}
    for target in graph.targets.values():
        hashes[target.name] = hash_spec(target.spec)
    return hashes


def _get_spec_hashes_path(working_dir):
    return os.path.join(working_dir, ".gwf", "spec-hashes.json")


def load_spec_hashes(working_dir, graph):
    path = _get_spec_hashes_path(working_dir)

    # If spec-hashes.json doesn't exist, it's because this is the first time
    # we're using spec hashes, so we'll initialize the dict with the current
    # spec hashes.
    if not os.path.exists(path):
        logging.debug("First run with spec hashes enabled - computing hashes")
        return hash_specs(graph)

    with open(path) as hashes_file:
        return json.load(hashes_file)


def dump_spec_hashes(working_dir, graph):
    hashes = hash_specs(graph)
    with open(_get_spec_hashes_path(working_dir), "w") as hashes_file:
        json.dump(hashes, hashes_file)


class TargetStatus(Enum):
    """Status of a target as computed by the Scheduler."""

    SHOULDRUN = 0  #: The target should run.
    SUBMITTED = 1  #: The target has been submitted, but is not currently running.
    RUNNING = 2  #: The target is currently running.
    COMPLETED = 3  #: The target has completed and should not run.


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
            self.protect = set()
        else:
            self.protect = set(protect)

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
        return (
            "{}(inputs={!r}, outputs={!r}, options={!r}, working_dir={!r}, "
            "spec={!r})"
        ).format(
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

    def protected(self):
        return set(_norm_paths(self.working_dir, _flatten(self.protect)))

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
        return "{}(name={!r})".format(self.__class__.__name__, self.name)

    def __str__(self):
        return self.name


class CircularDependencyError(GWFError):
    pass


class FileProvidedByMultipleTargetsError(GWFError):
    pass


@timer("Checked for circular dependencies in %.3fms", logger=logger)
def check_for_circular_dependencies(targets, dependencies):
    logger.debug("Checking for circular dependencies")

    fresh, started, done = 0, 1, 2

    nodes = targets.values()
    state = dict((n, fresh) for n in nodes)

    def visitor(node):
        state[node] = started
        for dep in dependencies[node]:
            if state[dep] == started:
                raise CircularDependencyError(
                    "Target {} depends on itself.".format(node)
                )
            elif state[dep] == fresh:
                visitor(dep)
        state[node] = done

    for node in nodes:
        if state[node] == fresh:
            visitor(node)


@attrs.define
class Graph:
    """Represents a dependency graph for a set of targets.

    The graph represents the targets present in a workflow, but also their
    dependencies and the files they provide.

    During construction of the graph the dependencies between targets are
    determined by looking at target inputs and outputs. If a target specifies a
    file as input, the file must either be provided by another target or
    already exist on disk. In case that the file is provided by another target,
    a dependency to that target will be added:

    :ivar dict dependencies:
        A dictionary mapping a target to a set of its dependencies.

    If the file is not provided by another target, the file is *unresolved*:

    :ivar set unresolved:
        A set containing file paths of all unresolved files.

    If the graph is constructed successfully, the following instance variables
    will be available:

    :ivar dict targets:
        A dictionary mapping target names to instances of :class:`gwf.Target`.
    :ivar dict provides:
        A dictionary mapping a file path to the target that provides that path.
    :ivar dict dependents:
        A dictionary mapping a target to a set of all targets which depend on
        the target.

    The graph can be manipulated in arbitrary, diabolic ways after it has been
    constructed. Checks are only performed at construction-time, thus
    introducing e.g. a circular dependency by manipulating *dependencies* will
    not raise an exception.
    """

    targets: dict
    provides: dict
    dependencies: defaultdict
    dependents: defaultdict
    unresolved: set

    @classmethod
    def from_targets(cls, targets):
        """Construct a dependency graph from a set of targets.

        When a graph is initialized it computes all dependency relations
        between targets, ensuring that the graph is semantically sane.
        Therefore, construction of the graph is an expensive operation which
        may raise a number of exceptions:

        :raises gwf.exceptions.FileProvidedByMultipleTargetsError:
            Raised if the same file is provided by multiple targets.

        :raises gwf.exceptions.CircularDependencyError:
            Raised if the graph contains a circular dependency.
        """
        provides = {}
        unresolved = set()
        dependencies = defaultdict(set)
        dependents = defaultdict(set)

        logger.debug("Building dependency graph from %d targets", len(targets))

        if isinstance(targets, dict):
            targets = targets.values()

        with timer("Built dependency graph in %.3fms", logger=logger):
            for target in targets:
                for path in target.flattened_outputs():
                    if path in provides:
                        msg = 'File "{}" provided by targets "{}" and "{}".'.format(
                            path, provides[path].name, target
                        )
                        raise FileProvidedByMultipleTargetsError(msg)
                    provides[path] = target

        for target in targets:
            for path in target.flattened_inputs():
                if path in provides:
                    dependencies[target].add(provides[path])
                else:
                    unresolved.add(path)

        for target, deps in dependencies.items():
            for dep in deps:
                dependents[dep].add(target)

        targets = {target.name: target for target in targets}
        check_for_circular_dependencies(targets, dependencies)

        return cls(
            targets=targets,
            provides=provides,
            dependencies=dependencies,
            dependents=dependents,
            unresolved=unresolved,
        )

    def endpoints(self):
        """Return a set of all targets that are not depended on by other targets."""
        return set(self.targets.values()) - set(self.dependents.keys())

    def dfs(self, root):
        """Return the depth-first traversal path through a graph from `root`."""
        visited = set()
        path = []

        def dfs_inner(node):
            if node in visited:
                return

            visited.add(node)
            for dep in self.dependencies[node]:
                dfs_inner(dep)
            path.append(node)

        dfs_inner(root)
        return path

    def __iter__(self):
        return iter(self.targets.values())

    def __len__(self):
        return len(self.targets)

    def __getitem__(self, target_name):
        return self.targets[target_name]

    def __contains__(self, target_name):
        return target_name in self.targets


@attrs.define
class CachedFilesystem:
    """A cached file system abstraction."""

    _cache: dict = attrs.field(factory=dict)

    def _lookup_file(self, path):
        if path not in self._cache:
            try:
                st = os.stat(path)
            except FileNotFoundError:
                self._cache[path] = None
            else:
                self._cache[path] = st.st_mtime
        return self._cache[path]

    def exists(self, path):
        return self._lookup_file(path) is not None

    def changed_at(self, path):
        st = self._lookup_file(path)
        if st is None:
            raise FileNotFoundError(path)
        return st
