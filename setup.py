# -*- coding: utf-8 -*-

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
            'run = gwf.plugins.run:run',
            'config = gwf.plugins.config:config',
            'status = gwf.plugins.status:status',
            'logs = gwf.plugins.logs:logs',
            'clean = gwf.plugins.clean:clean',
            'workers = gwf.plugins.workers:workers',
            'cancel = gwf.plugins.cancel:cancel'
        ]
    },

    test_suite='tests',
    setup_requires=['pytest-runner'],
    tests_require=['pytest',],
    install_requires=[
        'click>=6.6',
        'click-plugins>=0.2.2',
        'statusbar>=0.1.4',
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
        'Operating System :: OS Independent',
        'Topic :: Utilities',
        'Topic :: System :: Distributed Computing',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
)
