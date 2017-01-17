import os

from setuptools import find_packages, setup

import versioneer

# Utility function to read the README file.  Used for the
# long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put
# a raw string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='gwf',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'gwf = gwf.cli:main',
        ],
        'gwf.backends': [
            'slurm = gwf.backends.slurm:SlurmBackend',
            'local = gwf.backends.local:LocalBackend',
            'testing = gwf.backends.testing:TestingBackend',
        ],
        'gwf.plugins': [
            'run = gwf.plugins.run:RunCommand',
            'config = gwf.plugins.config:ConfigCommand',
            'status = gwf.plugins.status:StatusCommand',
            'logs = gwf.plugins.logs:LogsCommand',
            'clean = gwf.plugins.clean:CleanCommand',
            'workers = gwf.plugins.workers:WorkersPlugin'
        ]
    },

    test_suite='tests',
    install_requires=[
        'configargparse>=0.11',
        'statusbar>=0.1.2',
    ],

    # metadata for upload to PyPI
    author='Thomas Mailund, Dan Søndergaard',
    author_email='mailund@birc.au.dk, das@birc.au.dk',
    license='GPLv3',
    keywords='grid computing workflow',
    url='http://gwf.readthedocs.io/',
    description='A flexible, pragmatic workflow tool.',
    long_description=read('README.rst'),

    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Topic :: Utilities',
        'Topic :: System :: Distributed Computing',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
    ],
)
