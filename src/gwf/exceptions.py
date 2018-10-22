import click


class GWFError(click.ClickException):
    pass


class NameError(GWFError):
    pass


class TypeError(GWFError):
    pass


class WorkflowError(GWFError):
    pass


class ConfigurationError(GWFError):
    pass
