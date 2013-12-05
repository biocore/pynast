#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The PyNAST Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
from distutils.core import setup
import re

__version__ = "1.2.2"

# classes/classifiers code adapted from Celery and pyqi:
# https://github.com/celery/celery/blob/master/setup.py
# https://github.com/bipy/pyqi/blob/master/setup.py
#
# PyPI's list of classifiers can be found here:
# https://pypi.python.org/pypi?%3Aaction=list_classifiers
classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: User Interfaces 
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: OS Independent
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

# long_despcription should be all of the information from README.md, minus 
# the jenkins build status line
long_description = ''.join([line for line in open('README.md')
                            if not line.startswith("[![Build Status]")])

setup(name='pynast',
      version=__version__,
      description='The Python Nearest Alignment Space Termination tool',
      author="Greg Caporaso",
      author_email="gregcaporaso@gmail.com",
      maintainer="Greg Caporaso",
      maintainer_email="gregcaporaso@gmail.com",
      url='http://qiime.org/pynast',
      packages=['pynast','pynast/pycogent_backports'],
      scripts=['scripts/pynast'],
      long_description=long_description,
      install_requires=["numpy >= 1.5.1",
                        "cogent >= 1.5.3"],
      extras_require={'doc':"Sphinx >= 0.3"},
      classifiers=classifiers
)
