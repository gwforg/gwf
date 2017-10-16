from pkg_resources import iter_entry_points


from .base import (
    Backend,
    Status,
)


BACKENDS = {ep.name: ep.load() for ep in iter_entry_points('gwf.backends')}


def backend_from_config(config):
    return BACKENDS[config['backend']]