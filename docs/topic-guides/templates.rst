.. _templates:

=========
Templates
=========

Templates in *gwf* provide a simple mechanism for constructing a bunch of similar targets. For example, you may have
100 images that you want to transform in some way. We could write a *gwf* workflow to do this::

    from gwf import Workflow

    gwf = Workflow()

    photos = gwf.glob('photos/*.jpg')
    for index, path in enumerate(photos):
        gwf.target('TransformPhoto.{}'.format(index), inputs=[path], outputs=[path + '.new']) << """
        ./transform_photo {}
        """.format(path)

This will generate a target for each photo we've got which is all fine and dandy. However, we can make the code a lot
clearer using *templates*. In *gwf* a template is a function that returns a tuple containing four things:

1. The *input* files, corresponding to the *inputs* argument,
2. The *output* files, corresponding to the *outputs* argument,
3. a dictionary with options for the target that is to be generated, for example how many
   cores the template needs and which files it depends on,
4. a string which contains the specification of the target that is to be generated.

Let's rewrite our workflow from before with a template. First, we'll define the template function::

    def transform_photo(path):
        inputs = [path]
        outputs = [path + '.new']
        options = {}
        spec = """./transform_photo {}""".format(path)
        return inputs, outputs, options, spec

Next, we tell *gwf* to create a target using the :func:`~gwf.Workflow.target_from_template` method::

    from gwf import Workflow

    gwf = Workflow()

    photos = gwf.glob('photos/*.jpg')
    for index, path in enumerate(photos):
        gwf.target_from_template('TransformPhoto.{}'.format(index), transform_photo(path))

Templates are just Python functions, so you can do pretty much anything in a template function. For example, you can
create template functions that work across a wide range of systems. E.g. the template can determine the operating system
used and adapt the *spec* according to this.

Templates can also be put in separate files. We can put the ``transform_photo()`` template function into
``templates.py`` and then import it as we would import any other Python module::

    from gwf import Workflow
    from templates import transform_photo

    gwf = Workflow()

    photos = gwf.glob('photos/*.jpg')
    for index, path in enumerate(photos):
        gwf.target_from_template('TransformPhoto.{}'.format(index), transform_photo(path))

