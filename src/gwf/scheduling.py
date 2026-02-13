import logging
from functools import partial, lru_cache

from .backends.base import BackendStatus
from .core import Status

logger = logging.getLogger(__name__)


SUBMITTED_STATES = (
    Status.SUBMITTED,
    Status.RUNNING,
    Status.SHOULDRUN,
    Status.FAILED,
    Status.CANCELLED,
)


def schedule(
    endpoints,
    graph,
    fs,
    spec_hashes,
    status_func,
    submit_func,
    force=False,
    no_deps=False,
):
    @lru_cache(maxsize=None)
    def _non_temp_outputs(target):
        return frozenset(target.flattened_outputs()) - graph.temporary

    @lru_cache(maxsize=None)
    def _non_temp_inputs(target):
        return frozenset(target.flattened_inputs()) - graph.temporary

    @lru_cache(maxsize=None)
    def _oldest_output_mtime(target):
        outputs = _non_temp_outputs(target)
        missing = next((p for p in outputs if not fs.exists(p)), None)
        if missing is not None:
            return None, missing
        return min((fs.changed_at(p) for p in outputs), default=float("inf")), None

    @lru_cache(maxsize=None)
    def _youngest_input_mtime(target):
        inputs = _non_temp_inputs(target)
        first_missing = next((p for p in inputs if not fs.exists(p)), None)
        existing = [p for p in inputs if fs.exists(p)]
        youngest = max((fs.changed_at(p) for p in existing), default=float("-inf"))
        return youngest, first_missing

    @lru_cache(maxsize=None)
    def _dependents_are_up_to_date(target):
        if not (dependents := graph.dependents[target]):
            return False
        for dep in dependents:
            if _non_temp_outputs(dep):
                oldest_out, _ = _oldest_output_mtime(dep)
                if oldest_out is None:
                    return False
                youngest_in, _ = _youngest_input_mtime(dep)
                if youngest_in > oldest_out:
                    return False
            elif not _dependents_are_up_to_date(dep):
                return False
        return True

    @lru_cache(maxsize=None)
    def _should_run(target):
        if spec_hashes.has_changed(target) is not None:
            logger.debug("Target %s has a changed spec", target)
            return True

        # If I have no outputs, but I have inputs, I should probably only run if my input
        # changed, but I don't have any output files to compare with, so I'll just run
        # every time.
        if not target.outputs:
            logger.debug(
                "Target %s has no outputs and will always be scheduled", target
            )
            return True

        non_temp = _non_temp_outputs(target)
        temp_outputs = frozenset(target.flattened_outputs()) - non_temp

        # If I have non-temporary outputs, I should run if any of them are missing, or if any
        # existing input is newer than any existing output.
        if non_temp:
            oldest_out, missing_out = _oldest_output_mtime(target)
            if oldest_out is None:
                logger.debug("Target %s is missing output file %s", target, missing_out)
                return True

            youngest_in, missing_in = _youngest_input_mtime(target)
            if missing_in is not None:
                logger.debug("Target %s is missing input file %s", target, missing_in)
                return True

            if youngest_in > oldest_out:
                logger.debug("Target %s is not up-to-date", target)
                return True

            logger.debug("Target %s is up-to-date", target)
            return False

        # If I have only temporary outputs and any of them are missing, I should run if any downstream
        # targets are missing outputs or have newer inputs, otherwise I can consider myself up-to-date.
        missing_temp = next((p for p in temp_outputs if not fs.exists(p)), None)
        if missing_temp is not None:
            if _dependents_are_up_to_date(target):
                logger.debug(
                    "Target %s has missing temporary outputs but all downstream outputs exist and are up-to-date",
                    target,
                )
                return False

            logger.debug(
                "Target %s is missing temporary output file %s",
                target,
                missing_temp,
            )
            return True

        # If I have only temporary outputs and none of them are missing, I should run if any existing input
        # is newer than any existing output.
        oldest_temp = min(
            (fs.changed_at(p) for p in temp_outputs), default=float("inf")
        )
        if any(
            fs.exists(p) and fs.changed_at(p) > oldest_temp
            for p in target.flattened_inputs()
        ):
            logger.debug("Target %s (temp outputs) is not up-to-date", target)
            return True

        logger.debug("Target %s is up-to-date", target)
        return False

    def _schedule(target):
        submitted_deps = []

        if not no_deps:
            for dep in sorted(graph.dependencies[target], key=lambda t: t.name):
                status = _cached_schedule(dep)
                if status in SUBMITTED_STATES:
                    submitted_deps.append(dep)

        if force:
            logger.debug("Target %s is being forcibly submitted", target)
            submit_func(target, dependencies=submitted_deps)
            return Status.SHOULDRUN

        if status_func(target) == BackendStatus.SUBMITTED:
            logger.debug("Target %s is already submitted", target)
            return Status.SUBMITTED

        if status_func(target) == BackendStatus.RUNNING:
            logger.debug("Target %s is already running", target)
            return Status.RUNNING

        if status_func(target) == BackendStatus.FAILED:
            submit_func(target, dependencies=submitted_deps)
            return Status.FAILED

        if status_func(target) == BackendStatus.CANCELLED:
            submit_func(target, dependencies=submitted_deps)
            return Status.CANCELLED

        if submitted_deps:
            logger.debug(
                "Target %s will be submitted because of dependency %s",
                target,
                submitted_deps[0],
            )
            submit_func(target, dependencies=submitted_deps)
            return Status.SHOULDRUN

        if _should_run(target):
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

    # logger.info("Submitting target %s", target)

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

    if hasattr(backend, "get_tracked_id"):
        logger.info(
            "Submitted target %s (id: %s)",
            target,
            backend.get_tracked_id(target),
        )
    else:
        logger.info("Submitted target %s", target)


def submit_workflow(
    endpoints,
    graph,
    fs,
    spec_hashes,
    backend,
    dry_run=False,
    force=False,
    no_deps=False,
):
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
        force=force,
        no_deps=no_deps,
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
