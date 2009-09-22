#!/usr/bin/env python

import logging

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2009, the PyNAST project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Beta"

class NastLogger:
    __LABELS = [
        "candidate sequence ID", 
        "candidate nucleotide count",
        "errors", 
        "template ID", 
        "BLAST percent identity to template",
        "longest insertion relative to template",
        "candidate span aligned",
        "candidate nucleotide count post-NAST",
        "unaligned length",
        "count of single nucelotide 7mers or longer Nmers",
        "non-ACGT nucleotide count",
        "non-ACGT nucleotide percent",
        ]

    def __init__(self, filename=None):
        self.Filename = filename
        self.__logger = self.__init_logger()
        self.record(*self.__LABELS)

    def __init_logger(self):
        if self.Filename is not None:
            handler = logging.FileHandler(self.Filename, mode='w')
        else:
            class NullHandler(logging.Handler):
                def emit(self, record): pass
            handler = NullHandler()
        logger = logging.getLogger("PyNAST logger")
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def record(self, *args):
        log_entry = '\t'.join(map(str, args))
        self.__logger.info(log_entry)


