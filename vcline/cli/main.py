#!/usr/bin/env python
"""
Variant Calling Pipeline for Clinical Sequencing

Usage:
    vcline init [--debug|--info] [--yml=<path>]
    vcline download [--debug|--info] [--hg19] [<work_dir_path>]
    vcline run [--debug|--info] [--yml=<path>] [--cpus=<int>]
        [--ref-dir=<path>] [<work_dir_path>]
    vcline -h|--help
    vcline --version

Commands:
    init                Create a config YAML template
    download            Download hg38 reference FASTA
    run                 Run the analytical pipeline

Options:
    -h, --help          Print help and exit
    --version           Print version and exit
    --debug, --info     Execute a command with debug|info messages
    --cpus=<int>        Limit CPU cores used
    --yml=<path>        Specify a config YAML path [default: vcline.yml]
    --hg19              Use hg19 instead of hg38
    --ref-dir=<path>    Specify a reference directory path

Args:
    <work_dir_path>     Working directory path [default: .]
"""

import logging
import os
import shutil
from pathlib import Path

from docopt import docopt

from .. import __version__
from .pipeline import run_analytical_pipeline
from .util import download_and_merge_files, download_file, print_log, read_yml


def main():
    args = docopt(__doc__, version=__version__)
    _set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    if args['init']:
        _write_config_yml(path=args['--yml'])
    elif args['download']:
        _download_hg(
            hg_ver=('hg19' if args['--hg19'] else 'hg38'),
            work_dir_path=args['<work_dir_path>']
        )
    elif args['run']:
        run_analytical_pipeline(
            config_yml_path=args['--yml'], ref_dir_path=args['--ref-dir'],
            work_dir_path=args['<work_dir_path>'], n_cpu=args['--cpus']
        )


def _set_log_config(debug=None, info=None):
    if debug:
        lv = logging.DEBUG
    elif info:
        lv = logging.INFO
    else:
        lv = logging.WARNING
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S', level=lv
    )


def _write_config_yml(path):
    if Path(path).is_file():
        print_log('The file exists:\t{}'.format(path))
    else:
        print_log('Create a config YAML:\t{}'.format(path))
        shutil.copyfile(
            str(Path(__file__).parent.joinpath('../static/vcline.yml')),
            Path(path).resolve()
        )


def _download_hg(hg_ver='hg38', work_dir_path=None):
    print_log('Download genome FASTA:\t{}'.format(hg_ver))
    urls = read_yml(
        path=str(Path(__file__).parent.joinpath('../static/urls.yml'))
    )['ref_fa_gz'][hg_ver]
    assert all([u.endswith('.gz') for u in urls]), 'invalid gzip URLs'
    output_path = Path(work_dir_path or '.').joinpath(
        '.'.join([Path(Path(Path(u).name).stem).stem for u in urls])
        + '.fa.gz'
    )
    if len(urls) == 1:
        download_file(url=urls[0], output_path=output_path)
    else:
        download_and_merge_files(urls=urls, output_path=output_path)
