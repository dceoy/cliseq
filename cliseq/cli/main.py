#!/usr/bin/env python
"""
Analytical pipeline for clinical DNA sequencing data

Usage:
    cliseq init [--debug|--info] [--file=<yaml>]
    cliseq download [--debug|--info] [--hg19] [<dir>]
    cliseq run [--debug|--info] [--file=<yaml>] [--ref-dir=<dir>] [<dir>]
    cliseq -h|--help
    cliseq --version

Commands:
    init                Create a config YAML template
    download            Download hg38 reference FASTA
    run                 Run the analytical pipeline

Options:
    -h, --help          Print help and exit
    --version           Print version and exit
    --debug, --info     Execute a command with debug|info messages
    --file=<yaml>       Specify a config YAML path [default: cliseq.yml]
    --hg19              Use hg19 instead of hg38
    --ref-dir=<dir>     Specify a reference directory path

Args:
    <dir>               Working directory path [default: .]
"""

import logging
import os

from docopt import docopt

from . import __version__
from .pipeline import run_analytical_pipeline
from .util import download_hg, set_log_config, write_config_yml


def main():
    args = docopt(__doc__, version=__version__)
    set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    if args['init']:
        write_config_yml(path=args['--file'])
    elif args['download']:
        download_hg(
            ver=('hg19' if args['--hg19'] else 'hg38'), work_dir=args['<dir>']
        )
    elif args['run']:
        run_analytical_pipeline(
            config_yml_path=args['--file'], ref_dir_path=args['--ref-dir'],
            work_dir=args['<dir>']
        )
