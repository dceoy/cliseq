#!/usr/bin/env python
"""
Variant Calling Pipeline for Clinical Sequencing

Usage:
    vcline init [--debug|--info] [--yml=<path>]
    vcline run [--debug|--info] [--yml=<path>] [--cpus=<int>]
        [--workers=<int>] [--ref-dir=<path>] [<work_dir_path>]
    vcline -h|--help
    vcline --version

Commands:
    init                Create a config YAML template
    run                 Run the analytical pipeline

Options:
    -h, --help          Print help and exit
    --version           Print version and exit
    --debug, --info     Execute a command with debug|info messages
    --cpus=<int>        Limit CPU cores used
    --workers=<int>     Specify the maximum number of workers [default: 2]
    --yml=<path>        Specify a config YAML path [default: vcline.yml]
    --ref-dir=<path>    Specify a reference directory path

Args:
    <work_dir_path>     Working directory path [default: .]
"""

import logging
import os
import shutil
from pathlib import Path

import coloredlogs
from docopt import docopt

from .. import __version__
from .pipeline import run_analytical_pipeline
from .util import print_log


def main():
    args = docopt(__doc__, version=__version__)
    if args['--debug']:
        log_level = 'DEBUG'
    elif args['--info']:
        log_level = 'INFO'
    else:
        log_level = 'WARNING'
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S', level=log_level
    )
    coloredlogs.install(level=log_level)
    logger = logging.getLogger(__name__)
    logger.debug(f'args:{os.linesep}{args}')
    if args['init']:
        _write_config_yml(path=args['--yml'])
    elif args['run']:
        run_analytical_pipeline(
            config_yml_path=args['--yml'], ref_dir_path=args['--ref-dir'],
            work_dir_path=args['<work_dir_path>'], max_n_cpu=args['--cpus'],
            max_n_worker=args['--workers'], log_level=log_level
        )


def _write_config_yml(path):
    if Path(path).is_file():
        print_log(f'The file exists:\t{path}')
    else:
        print_log(f'Create a config YAML:\t{path}')
        shutil.copyfile(
            str(Path(__file__).parent.joinpath('../static/vcline.yml')),
            Path(path).resolve()
        )
