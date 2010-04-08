#!/usr/bin/env python

import logging

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2010, The PyNAST Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Release"

class NastLogger:
    __LABELS = [
        "candidate sequence ID", 
        "candidate nucleotide count",
        "errors", 
        "template ID", 
        "BLAST percent identity to template",
        "candidate nucleotide count post-NAST",
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


