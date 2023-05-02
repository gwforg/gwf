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
import click

from .exceptions import GWFError
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


class InvalidPathError(GWFError):
    pass


def _check_path(path):
    if not path:
        raise InvalidPathError("Path is empty")

    result = _has_nonprintable_char(path)
    if result is not None:
        clean_path, char, pos = result
        raise InvalidPathError(
            f"Path {clean_path} contains a non-printable character '{char}' at {pos}. "
            f"This is always unintentional and can cause strange behaviour."
        )


def _norm_path(working_dir, path):
    path = fspath(path)
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(working_dir, path))


def _norm_paths(working_dir, paths):
    return [_norm_path(working_dir, path) for path in paths]


def hash_spec(spec):
    return hashlib.sha1(spec.encode("utf-8")).hexdigest()


@attrs.define
class NoopSpecHashes:
    def has_changed(self, target):  # pragma: no cover
        return None

    def update(self, target):  # pragma: no cover
        pass

    def invalidate(self, target):  # pragma: no cover
        pass

    def close(self):  # pragma: no cover
        pass

    def __enter__(self):  # pragma: no cover
        return self

    def __exit__(self, *exc):  # pragma: no cover
        self.close()


@attrs.define
class FileSpecHashes:
    path: str = attrs.field()
    hashes: dict = attrs.field(factory=dict, init=False)

    def __attrs_post_init__(self):
        try:
            with open(self.path) as hashes_file:
                self.hashes = json.load(hashes_file)
        except FileNotFoundError:
            logger.debug("First run with spec hashes enabled")

    def has_changed(self, target):
        spec_hash = hash_spec(target.spec)
        saved_hash = self.hashes.get(target.name)
        if saved_hash is None:
            logger.debug("No spec hash for %s exists", target)
            return spec_hash
        if spec_hash != saved_hash:
            return spec_hash
        return None

    def update(self, target):
        self.hashes[target.name] = hash_spec(target.spec)

    def invalidate(self, target):
        try:
            del self.hashes[target.name]
        except KeyError:
            pass

    def close(self):
        with open(self.path, "w") as hashes_file:
            json.dump(self.hashes, hashes_file)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def get_spec_hashes(*, working_dir, config):
    if config.get("use_spec_hashes"):
        return FileSpecHashes(os.path.join(working_dir, ".gwf", "spec-hashes.json"))
    else:
        return NoopSpecHashes()


class Status(Enum):
    """BackendStatus of a target as computed by the Scheduler."""

    SHOULDRUN = 0  #: The target should run.
    SUBMITTED = 1  #: The target has been submitted, but is not currently running.
    RUNNING = 2  #: The target is currently running.
    COMPLETED = 3  #: The target has completed and should not run.
    FAiLED = 4  #: The target failed and should run again.


@attrs.define(eq=False)
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

    inputs: list = attrs.field()
    outputs: list = attrs.field()
    options: dict = attrs.field()
    working_dir: str = attrs.field(default=".")
    protect: set = attrs.field(factory=set, converter=set)
    spec: str = attrs.field(default="")


@attrs.define(eq=False)
class Target:
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

    name: str = attrs.field()
    inputs: list = attrs.field()
    outputs: list = attrs.field()
    options: dict = attrs.field()
    working_dir: str = attrs.field(default=".")
    protect: set = attrs.field(factory=set, converter=set)
    spec: str = attrs.field(default="")
    order: int = attrs.field(init=False)

    _creation_order = 0

    @order.default  # type: ignore
    def _set_target_order(self):
        order = Target._creation_order
        Target._creation_order += 1
        return order

    @name.validator  # type: ignore
    def _validate_name(self, attribute, value):
        if not is_valid_name(self.name):
            raise GWFError(f"Target defined with invalid name: {value}")

    @inputs.validator
    @outputs.validator
    def _validate_inputs(self, attribute, value):
        for path in value:
            _check_path(path)

    @working_dir.validator
    def _validate_working_dir(self, attribute, value):
        _check_path(value)

    def __lshift__(self, spec):
        self.spec = spec
        return self

    def __str__(self):
        return self.name

    def flattened_inputs(self):
        return _norm_paths(self.working_dir, _flatten(self.inputs))

    def flattened_outputs(self):
        return _norm_paths(self.working_dir, _flatten(self.outputs))

    def protected(self):
        return set(_norm_paths(self.working_dir, _flatten(self.protect)))


class CircularDependencyError(GWFError):
    pass


class FileProvidedByMultipleTargetsError(GWFError):
    pass


class UnresolvedInputError(GWFError):
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

    targets: dict = attrs.field()
    provides: dict = attrs.field()
    dependencies: defaultdict = attrs.field()
    dependents: defaultdict = attrs.field()
    unresolved: set = attrs.field()

    @classmethod
    def from_targets(cls, targets, fs):
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

        # Check whether all input files actually exist or are being provided
        # by another target. If not, it's an error.
        for target in targets:
            for path in target.flattened_inputs():
                if path in unresolved and not fs.exists(path):
                    msg = (
                        'File "{}" is required by "{}", but does not exist and is not '
                        "provided by any target in the workflow."
                    ).format(path, target)
                    raise UnresolvedInputError(msg)

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


@attrs.frozen
class Context:
    working_dir = attrs.field()
    config = attrs.field()
    backend = attrs.field()
    workflow_file = attrs.field()
    workflow_obj = attrs.field()

    @property
    def config_dir(self):
        return os.path.join(self.working_dir, ".gwf")

    @property
    def logs_dir(self):
        return os.path.join(self.config_dir, "logs")


pass_context = click.make_pass_decorator(Context)
