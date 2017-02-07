import os.path

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from gwf.core import PreparedWorkflow, Target, Workflow  # noqa: E402

USER_CONFIG_FILE = os.path.expanduser('~/.gwf.conf')
LOCAL_CONFIG_FILE = '.gwf.conf'


class template(object):
    """A template with string substitution.

    .. deprecated:: 1.0

        Use function templates as described in :ref:`function_templates`.
    """

    def __init__(self, **options):
        self.spec = None
        self.options = options

    def __lshift__(self, spec):
        self.spec = spec
        return self

    def __call__(self, **substitutions):
        def substitute(s):
            if isinstance(s, str):
                return s.format(**substitutions)
            elif hasattr(s, '__iter__'):
                return [substitute(x) for x in s]

        formatted_options = [(key, substitute(val))
                             for key, val in self.options.items()]
        options = dict(formatted_options)
        return options, self.spec.format(**substitutions)

    def __repr__(self):
        return '{}(options={!r}, spec={!r})'.format(
            self.__class__.__name__, self.options, self.spec)


__all__ = ('Target', 'Workflow', 'PreparedWorkflow', 'template',)
