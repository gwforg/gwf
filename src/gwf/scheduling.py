import logging
from functools import partial

from .backends.base import BackendStatus
from .core import Status

logger = logging.getLogger(__name__)


SUBMITTED_STATES = (
    Status.SUBMITTED,
    Status.RUNNING,
    Status.SHOULDRUN,
)


def should_run(target, fs, spec_hashes):
    new_hash = spec_hashes.has_changed(target)
    if new_hash is not None:
        logger.debug("Target %s has a changed spec", target)
        return True

    for path in target.flattened_outputs():
        if not fs.exists(path):
            logger.debug("Target %s is missing output file %s", target, path)
            return True

    if not target.inputs:
        logger.debug("Target %s has no inputs", target)
        return True

    youngest_in_ts, youngest_in_path = max(
        (fs.changed_at(path), path)
        for path in target.flattened_inputs()
        if fs.exists(path)
    )

    # If I have no outputs, but I have inputs, I should probably only run if my input
    # changed, but I don't have any output files to compare with, so I'll just run
    # every time.
    if not target.outputs:
        logger.debug("Target %s has no outputs", target)
        return True

    oldest_out_ts, oldest_out_path = min(
        (fs.changed_at(path), path)
        for path in target.flattened_outputs()
        if fs.exists(path)
    )

    if youngest_in_ts > oldest_out_ts:
        logger.debug("Target %s is not up-to-date", target)
        return True

    logger.debug("Target %s is up-to-date", target)
    return False


def schedule(endpoints, graph, fs, spec_hashes, status_func, submit_func):
    def _schedule(target):
        submitted_deps = []
        for dep in graph.dependencies[target]:
            status = _cached_schedule(dep)
            if status in SUBMITTED_STATES:
                submitted_deps.append(dep)

        if status_func(target) == BackendStatus.SUBMITTED:
            logger.debug("Target %s is already submitted", target)
            return Status.SUBMITTED

        if status_func(target) == BackendStatus.RUNNING:
            logger.debug("Target %s is already running", target)
            return Status.RUNNING

        if status_func(target) == BackendStatus.FAILED:
            submit_func(target, dependencies=submitted_deps)
            return Status.FAiLED

        if submitted_deps:
            logger.debug(
                "Target %s will be submitted because of dependency %s",
                target,
                submitted_deps[0],
            )
            submit_func(target, dependencies=submitted_deps)
            return Status.SHOULDRUN

        if should_run(target, fs, spec_hashes):
            submit_func(target, dependencies=submitted_deps)
            return Status.SHOULDRUN

        return Status.COMPLETED

    cache = {}

    def _cached_schedule(target):
        if target not in cache:
            cache[target] = _schedule(target)
        return cache[target]

    for target in sorted(endpoints, key=lambda t: t.name):
        _cached_schedule(target)

    return cache


def _submit_dryrun(target, dependencies, backend, spec_hashes):
    logger.info("Would submit %s", target)


def _submit_noop(target, dependencies, backend, spec_hashes):
    pass


def submit_backend(target, dependencies, backend, spec_hashes):
    """Prepare and submit `target` with `dependencies` to `backend`.

    Will prepare the target for submission by injecting option defaults from the
    backend, check for unsupported options, and removing options with a `None`
    value.

    This is the primary way to submit a target. Do not call :func:`submit`
    directly on the backend, unless you want to manually deal with with
    injection of option defaults.
    """

    logger.info("Submitting target %s", target)

    new_options = {}
    if hasattr(backend, "target_defaults"):
        new_options = dict(backend.target_defaults)
    new_options.update(target.options)

    for option_name, option_value in list(new_options.items()):
        if option_name not in backend.target_defaults.keys():
            logger.warning(
                "Option '%s' used in '%s' is not supported by backend. Ignored.",
                option_name,
                target.name,
            )
            del new_options[option_name]
        elif option_value is None:
            del new_options[option_name]
    target.options = new_options

    backend.submit(target, dependencies)
    spec_hashes.update(target)


def submit_workflow(endpoints, graph, fs, spec_hashes, backend, dry_run=False):
    """Submit a workflow to a backend."""
    submit_func = partial(
        _submit_dryrun if dry_run else submit_backend,
        backend=backend,
        spec_hashes=spec_hashes,
    )
    schedule(
        endpoints,
        graph,
        fs,
        spec_hashes,
        status_func=backend.status,
        submit_func=submit_func,
    )


def get_status_map(graph, fs, spec_hashes, backend, endpoints=None):
    """Get the status of each targets in the graph."""
    submit_func = partial(_submit_noop, backend=backend, spec_hashes=spec_hashes)
    return schedule(
        endpoints if endpoints is not None else graph.endpoints(),
        graph,
        fs,
        spec_hashes,
        status_func=backend.status,
        submit_func=submit_func,
    )
