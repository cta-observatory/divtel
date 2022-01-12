#!/usr/bin/env python
# Licensed under an MIT license - see LICENSE

from setuptools import setup, find_packages
import os
import re


def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(project + '/__init__.py').read())
    return result.group(1)

long_description = open('README.md').read()

setup(name='divtel',
      version=get_property('__version__', 'divtel'),
      description="Divergent pointing mode for Imaging Atmospheric Cherenkov Telescopes arrays",
      packages=find_packages(),
      install_requires=['astropy',
                        'numpy'
                        ],
      tests_require=['pytest',
                     'pytest-ordering',
                     ],
      package_data={},
      author='T. Vuillaume, A. Donini, T. Gasparetto',
      author_email='thomas.vuillaume[at]lapp.in2p3.fr',
      license='MIT',
      url='https://github.com/cta-observatory/divtel',
      long_description=long_description,
      )
