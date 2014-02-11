#!/usr/bin/env python

from distutils.core import setup

setup(
    name='primer-lint',
    version='1.0.0',
    author='Luke Shillabeer',
    author_email='lshillabeer@gmail.com',
    packages=['primerlint'],
    scripts=[
        'primerlint/primerlint.py',
        'primerlint/importfiles.py',
        'primerlint/outputfiles.py',
        'primerlint/analysismodules.py'
        ],
    entry_points={
        'console_scripts': ['primer-lint = primerlint:main']
        },
    description=(
        'Primer-Lint: A tool for analysing primers for the Hi-Plex \
         targeted, multiplexed DNA sequencing strategy.'
         ),
    install_requires=['biopython >= 1.62'],
)