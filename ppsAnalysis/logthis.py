#!/usr/bin/env python#VERSION#

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import logging
import sys

class logit(object):

    def __init__(self, format=logging.Formatter("%(asctime)s — %(name)s — %(levelname)s — %(message)s"), log_f="", log_level=""):
        self._format = format
        self._logfile = log_f
        self._log_level = log_level

    def get_console_handler(self):
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(self._format)
        return console_handler

    def get_file_handler(self):
        file_handler = logging.FileHandler(self._logfile)
        file_handler.setFormatter(self._format)
        return file_handler

    def get_logger(self, logger_name):
        logger = logging.getLogger(logger_name)
        logger.setLevel(self._log_level)  # better to have too much log than not enough
        logger.addHandler(self.get_console_handler())
        logger.addHandler(self.get_file_handler())
        # with this pattern, it's rarely necessary to propagate the error up to parent
        logger.propagate = False
        return logger