#!/usr/bin/env python
"""
Variant Calling Pipeline for Clinical Sequencing

Usage:
    vcline init [--debug|--info] [--yml=<path>]
    vcline run [--debug|--info] [--yml=<path>] [--cpus=<int>]
        [--workers=<int>] [--split-intervals] [--ref-dir=<path>]
        [<work_dir_path>]
    vcline download-gnomad-af-vcf [--debug|--info] [--cpus=<int>]
        [--vcf-bgz=<url>] [<work_dir_path>]
    vcline -h|--help
    vcline --version

Commands:
    init                Create a config YAML template
    run                 Run the analytical pipeline
    download-gnomad-af-vcf
                        Download a large VCF and extract only AF from INFO

Options:
    -h, --help          Print help and exit
    --version           Print version and exit
    --debug, --info     Execute a command with debug|info messages
    --yml=<path>        Specify a config YAML path [default: vcline.yml]
    --cpus=<int>        Limit CPU cores used
    --workers=<int>     Specify the maximum number of workers [default: 2]
    --split-intervals   Split evaluation intervals
    --ref-dir=<path>    Specify a reference directory path
    --vcf-bgz=<url>     Specify the URL of a bgzipped VCF

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
from .downloader import download_gnomad_vcf_and_extract_af
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
            max_n_worker=args['--workers'],
            split_intervals=args['--split-intervals'], log_level=log_level
        )
    elif args['download-gnomad-af-vcf']:
        download_gnomad_vcf_and_extract_af(
            vcf_bgz_url=args['--vcf-bgz'],
            work_dir_path=args['<work_dir_path>'], max_n_cpu=args['--cpus']
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
