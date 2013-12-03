#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The PyNAST Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
from tempfile import NamedTemporaryFile
from os import remove
from cogent import LoadSeqs, DNA
from cogent.util.unit_test import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from pynast.logger import NastLogger
from pynast.util import get_pynast_temp_dir

class NastLoggerTests(TestCase):
    """Tests of the PyNAST logging class"""

    def setUp(self):
        # Note that delete = False here because we don't want these to 
        # be deleted when they are closed (since we need to pass
        # the filepaths around after we write and close them). The files
        # are deleted explicitly at the end of the test.
        self.file = NamedTemporaryFile(prefix='NastLoggerTest',
                                       suffix='.log',
                                       dir=get_pynast_temp_dir(),
                                       delete=False)
        self.file.close()
        self.filename = self.file.name

    def tearDown(self):
        try:
            remove(self.filename)
        except OSError:
            pass

    def test_init(self):
        """NastLogger.__init__ should store log filename in Filename attribute"""
        null_logger = NastLogger()
        self.assertEqual(null_logger.Filename, None)

        file_logger = NastLogger(self.filename)
        self.assertEqual(file_logger.Filename, self.filename)

    def test_header(self):
        """NastLogger.__init__ should write correct header to log file"""
        logger = NastLogger(self.filename)

        f = open(self.filename, 'r')
        header = f.readline()
        f.close()

        exp_header = (
            'candidate sequence ID\tcandidate nucleotide count\terrors\t'
            'template ID\tBLAST percent identity to template\t'
            'candidate nucleotide count post-NAST\n'
            )
        self.assertEqual(header, exp_header)

    def test_record(self):
        """NastLogger.__init__ should record tab-separated values to log file"""

        logger = NastLogger(self.filename)
        logger.record('hello', 'world')

        f = open(self.filename, 'r')
        obs_header = f.readline()
        obs_message = f.readline()
        f.close()

        self.assertEqual(obs_message, 'hello\tworld\n')
        

if __name__ == "__main__":
    main()

