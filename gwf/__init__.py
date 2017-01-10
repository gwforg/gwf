from gwf.core import Event, PreparedWorkflow, Target, Workflow

USER_CONFIG_FILE = '~/.gwfrc'
LOCAL_CONFIG_FILE = '.gwfrc'


class template(object):
    """A template with string substitution."""

    def __init__(self, **options):
        self.spec = None
        self.options = options

    def __lshift__(self, spec):
        self.spec = spec
        return self

    def __call__(self, **substitutions):
        def substitute(s):
            if type(s) == str:
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
