#!/usr/bin/env python

from distutils.core import setup

setup(
    name='primer-lint',
    version='1.0.0',
    author='Luke Shillabeer',
    author_email='lshillabeer@gmail.com',
    packages=['primerlint'],
    scripts=['Primer_design/Primer_design.py'],
    entry_points={
        'console_scripts': ['primer-lint = primerlint:main']
    },
    url='https://github.com/Kiwwa/primer-lint',
    licence='LICENSE.txt',
    description=(
        'Hiplex-Primer: A tool for generating primers for the Hi-Plex \
         targeted, multiplexed DNA sequencing strategy.'),
    install_requires=['biopython >= 1.62'],
)