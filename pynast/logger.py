#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The PyNAST Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import logging

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


