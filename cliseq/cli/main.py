#!/usr/bin/env python
"""
Analytical pipeline for clinical DNA sequencing data

Usage:
    cliseq init [--debug|--info] [--yml=<path>]
    cliseq download [--debug|--info] [--hg19] [<work_dir_path>]
    cliseq run [--debug|--info] [--yml=<path>] [--ref-dir=<path>]
        [<work_dir_path>]
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
    --yml=<path>        Specify a config YAML path [default: cliseq.yml]
    --hg19              Use hg19 instead of hg38
    --ref-dir=<path>    Specify a reference directory path

Args:
    <work_dir_path>     Working directory path [default: .]
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
        write_config_yml(path=args['--yml'])
    elif args['download']:
        download_hg(
            hg_ver=('hg19' if args['--hg19'] else 'hg38'),
            work_dir_path=args['<work_dir_path>']
        )
    elif args['run']:
        run_analytical_pipeline(
            config_yml_path=args['--yml'], ref_dir_path=args['--ref-dir'],
            work_dir_path=args['<work_dir_path>']
        )
