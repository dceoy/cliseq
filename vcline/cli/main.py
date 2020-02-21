#!/usr/bin/env python
"""
Variant Calling Pipeline for Clinical Sequencing

Usage:
    vcline init [--debug|--info] [--yml=<path>]
    vcline run [--debug|--info] [--yml=<path>] [--cpus=<int>]
        [--workers=<int>] [--split-intervals] [--ref-dir=<path>]
        [<work_dir_path>]
    vcline download-funcotator-data [--debug|--info] [--cpus=<int>]
        [<work_dir_path>]
    vcline write-af-only-vcf [--debug|--info] [--cpus=<int>]
        [--src-path=<path>|--src-url=<url>|--gnomad-url] [<work_dir_path>]
    vcline -h|--help
    vcline --version

Commands:
    init                Create a config YAML template
    run                 Run the analytical pipeline
    download-funcotator-data
                        Download Funcotator data sources
    write-af-only-vcf   Extract and write only AF from VCF INFO

Options:
    -h, --help          Print help and exit
    --version           Print version and exit
    --debug, --info     Execute a command with debug|info messages
    --yml=<path>        Specify a config YAML path [default: vcline.yml]
    --cpus=<int>        Limit CPU cores used
    --workers=<int>     Specify the maximum number of workers [default: 2]
    --split-intervals   Split evaluation intervals
    --ref-dir=<path>    Specify a reference directory path
    --src-url=<url>     Specify a source URL
    --src-path=<path>   Specify a source PATH
    --gnomad-url        Specify the URL of gnomAD VCF (default)

Args:
    <work_dir_path>     Working directory path [default: .]
"""

import logging
import os
import shutil
from pathlib import Path

import coloredlogs
import luigi
from docopt import docopt
from psutil import cpu_count

from .. import __version__
from ..task.resource import DownloadFuncotatorDataSources, WriteAfOnlyVCF
from .pipeline import run_analytical_pipeline
from .util import fetch_executable, load_default_url_dict, print_log


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
    else:
        n_cpu = int(args['--cpus'] or cpu_count())
        if args['download-funcotator-data']:
            luigi.build(
                [
                    DownloadFuncotatorDataSources(
                        dest_dir_path=str(
                            Path(args['<work_dir_path>'] or '.').resolve()
                        ),
                        n_cpu=n_cpu, gatk=fetch_executable('gatk')
                    )
                ],
                local_scheduler=True, log_level=log_level
            )
        elif args['write-af-only-vcf']:
            luigi.build(
                [
                    WriteAfOnlyVCF(
                        src_path=(
                            str(Path(args['--src-path']).resolve())
                            if args['--src-path'] else ''
                        ),
                        src_url=(
                            '' if args['--src-path'] else (
                                args['--src-url']
                                or load_default_url_dict()['gnomad_vcf']
                            )
                        ),
                        dest_dir_path=str(
                            Path(args['<work_dir_path>'] or '.').resolve()
                        ),
                        n_cpu=n_cpu,
                        **{c: fetch_executable(c) for c in ['curl', 'bgzip']}
                    )
                ],
                local_scheduler=True, log_level=log_level
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
