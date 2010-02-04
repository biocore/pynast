#!/usr/bin/env python
# File created on 04 Feb 2010
from __future__ import division
from distutils.core import setup

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The PyNAST project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"
 
long_description = """The Python Nearest Alignment Search Tool
http://pynast.sourceforge.net

PyNAST: a flexible tool for aligning sequences to a template alignment. 
J. Gregory Caporaso, Kyle Bittinger, Frederic D. Bushman, Todd Z. DeSantis, Gary L. Andersen, and Rob Knight. 
January 15, 2010, DOI 10.1093/bioinformatics/btp636. Bioinformatics 26: 266-267.

"""
 
setup(name='PyNAST',
      version=__version__,
      description='The Python Nearest Alignment Search Tool',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://pynast.sourceforge.net',
      packages=['pynast'],
      scripts=['scripts/pynast'],
      long_description=long_description,
    )
