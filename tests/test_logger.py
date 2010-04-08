#!/usr/bin/env python

from __future__ import division
from os import remove
from cogent import LoadSeqs, DNA
from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from pynast.logger import NastLogger

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2010, The PyNAST Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Release"

class NastLoggerTests(TestCase):
    """Tests of the PyNAST logging class"""

    def setUp(self):
        self.filename = get_tmp_filename(
            prefix='NastLoggerTest',
            suffix='.log',
            )

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

        file = open(self.filename, 'r')
        header = file.readline()
        file.close()

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

        file = open(self.filename, 'r')
        obs_header = file.readline()
        obs_message = file.readline()
        file.close()

        self.assertEqual(obs_message, 'hello\tworld\n')
        

if __name__ == "__main__":
    main()

