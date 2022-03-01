#!/usr/bin/env python
# Licensed under an MIT license - see LICENSE

from setuptools import setup, find_packages
import os
import re


## read __version__
with open('divtel/version.py') as f:
    exec(f.read())

long_description = open('README.md').read()

setup(name='divtel',
      version=__version__,
      description="Divergent pointing mode for Imaging Atmospheric Cherenkov Telescopes arrays",
      packages=find_packages(),
      install_requires=['astropy',
                        'numpy',
                        'matplotlib',
                        ],
      extras_require={'tests': ['pytest', 'pytest-ordering']},
      package_data={},
      author='T. Vuillaume, A. Donini, T. Gasparetto',
      author_email='thomas.vuillaume[at]lapp.in2p3.fr',
      license='MIT',
      url='https://github.com/cta-observatory/divtel',
      long_description=long_description,
      )
