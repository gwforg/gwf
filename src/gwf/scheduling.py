import logging
from functools import partial

from .backends.base import Status
from .core import TargetStatus

logger = logging.getLogger(__name__)


SUBMITTED_STATES = (
    TargetStatus.SUBMITTED,
    TargetStatus.RUNNING,
    TargetStatus.SHOULDRUN,
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

        if status_func(target) == Status.SUBMITTED:
            logger.debug("Target %s is already submitted", target)
            return TargetStatus.SUBMITTED

        if status_func(target) == Status.RUNNING:
            logger.debug("Target %s is already running", target)
            return TargetStatus.RUNNING

        if status_func(target) == Status.FAILED:
            submit_func(target, dependencies=submitted_deps)
            return TargetStatus.FAiLED

        if submitted_deps:
            logger.debug(
                "Target %s will be submitted because of dependency %s",
                target,
                submitted_deps[0],
            )
            submit_func(target, dependencies=submitted_deps)
            return TargetStatus.SHOULDRUN

        if should_run(target, fs, spec_hashes):
            submit_func(target, dependencies=submitted_deps)
            return TargetStatus.SHOULDRUN

        return TargetStatus.COMPLETED

    cache = {}

    def _cached_schedule(target):
        if target not in cache:
            cache[target] = _schedule(target)
        return cache[target]

    for target in sorted(endpoints, key=lambda t: t.name):
        _cached_schedule(target)

    return cache


def _submit_dryrun(target, dependencies):
    logger.info("Would submit %s", target)


def submit_backend(target, dependencies, backend, spec_hashes):
    logger.info("Submitting target %s", target)
    backend.submit(target, dependencies)
    spec_hashes.update(target)


def _submit_noop(target, dependencies):
    pass


def submit_workflow(endpoints, graph, fs, spec_hashes, backend, dry_run=False):
    """Submit a workflow to a backend."""
    submit_func = (
        _submit_dryrun
        if dry_run
        else partial(submit_backend, backend=backend, spec_hashes=spec_hashes)
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
    return schedule(
        endpoints if endpoints is not None else graph.endpoints(),
        graph,
        fs,
        spec_hashes,
        status_func=backend.status,
        submit_func=_submit_noop,
    )
