#!/usr/bin/env python
# File created on 04 Feb 2010
from __future__ import division
from distutils.core import setup
import re

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The PyNAST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 
long_description = """The Python Nearest Alignment Space Termination tool
http://qiime.org/pynast

PyNAST: a flexible tool for aligning sequences to a template alignment. 
J. Gregory Caporaso, Kyle Bittinger, Frederic D. Bushman, Todd Z. DeSantis, Gary L. Andersen, and Rob Knight. 
January 15, 2010, DOI 10.1093/bioinformatics/btp636. Bioinformatics 26: 266-267.

"""
try:
    import cogent
except ImportError:
    print "PyCogent not installed but required. (Is it installed? Is it in the current users $PYTHONPATH or site-packages?) See http://pycogent.sourceforge.net."
    exit(1)

pycogent_version = tuple([int(v) \
        for v in re.split("[^\d]", cogent.__version__) if v.isdigit()])

if pycogent_version < (1,5,3):
    print "PyCogent >= 1.5.3 required, but %s is installed." % cogent.__version__
    exit(1)

setup(name='PyNAST',
      version=__version__,
      description='The Python Nearest Alignment Space Termination tool',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://qiime.org/pynast',
      packages=['pynast'],
      scripts=['scripts/pynast'],
      long_description=long_description
)
