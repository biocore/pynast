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

try:
    import cogent
except ImportError:
    print "PyCogent not installed but required. (Is it installed? Is it in the current user's $PYTHONPATH or site-packages?) See http://www.pycogent.org"
    exit(1)

pycogent_version = tuple([int(v) \
        for v in re.split("[^\d]", cogent.__version__) if v.isdigit()])

if pycogent_version < (1,5,3):
    print "PyCogent >= 1.5.3 required, but %s is installed." % cogent.__version__
    exit(1)

# long_despcription should be all of the information from README.md, minus 
# the jenkins build status line
long_description = ''.join([line for line in open('README.md')
                            if not line.startswith("[![Build Status]")])

setup(name='pynast',
      version="1.2.1-dev",
      description='The Python Nearest Alignment Space Termination tool',
      author="Greg Caporaso",
      author_email="gregcaporaso@gmail.com",
      maintainer="Greg Caporaso",
      maintainer_email="gregcaporaso@gmail.com",
      url='http://qiime.org/pynast',
      packages=['pynast'],
      scripts=['scripts/pynast'],
      long_description=long_description,
      install_requires=["cogent >= 1.5.3"],
      extras_require={'test':["nose >= 0.10.1",
                              "tox >= 1.6.1"],
                      'doc':"Sphinx >= 0.3"
                     },
      classifiers=classifiers
)
