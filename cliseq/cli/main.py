#!/usr/bin/env python
"""
Analytical pipeline for clinical DNA sequencing data

Usage:
    cliseq init [--debug|--info] [--file=<yaml>]
    cliseq run [--debug|--info] [--file=<yaml>] [<dir>]
    cliseq -h|--help
    cliseq --version

Options:
    -h, --help          Print help and exit
    --version           Print version and exit
    --debug, --info     Execute a command with debug|info messages
    --file=<yaml>       Specify a config YAML path [default: cliseq.yml]

Args:
    <dir>               Working directory [default: .]
"""

import logging
import os

from docopt import docopt

from . import __version__
from .util import set_log_config


def main():
    args = docopt(__doc__, version=__version__)
    set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    print(args)
