"""Module providing an alternative syntax for specifying targets."""


class target:
    """Class used to specify targets via operators instead of a
    method call to the workflow."""
    def __init__(self, name):
        self.name = name

    def __rfloordiv__(self, workflow):
        return workflow.target(self.name, inputs=[], outputs=[])


class ExplicitOption:
    """Abstract class for specifying options as explicit calls."""

    def update_target(self, target):
        """Overload this method to define a new option."""
        pass

    def __rfloordiv__(self, target):
        self.update_target(target)
        return target


class inputs(ExplicitOption):
    def __init__(self, *inputs):
        self.inputs = inputs

    def update_target(self, target):
        target.inputs = self.inputs


class outputs(ExplicitOption):
    def __init__(self, *outputs):
        self.outputs = outputs

    def update_target(self, target):
        target.outputs = self.outputs
