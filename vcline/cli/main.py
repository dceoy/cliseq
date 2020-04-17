#!/usr/bin/env python
"""
Variant Calling Pipeline for Clinical Sequencing

Usage:
    vcline init [--debug|--info] [--yml=<path>]
    vcline run [--debug|--info] [--yml=<path>] [--cpus=<int>]
        [--workers=<int>] [--split-intervals] [--skip-cleaning]
        [--ref-dir=<path>] [--dest-dir=<path>]
    vcline download-resources [--debug|--info] [--cpus=<int>]
        [--without-gnomad] [--dest-dir=<path>]
    vcline download-funcotator-data [--debug|--info] [--cpus=<int>]
        [--dest-dir=<path>]
    vcline download-gnomad-af-only-vcf [--debug|--info] [--cpus=<int>]
        [--dest-dir=<path>]
    vcline write-af-only-vcf [--debug|--info] [--cpus=<int>]
        [--src-path=<path>|--src-url=<url>] [--dest-dir=<path>]
    vcline create-interval-list [--debug|--info] [--cpus=<int>]
        [--dest-dir=<path>] <fa_path> <bed_path>
    vcline funcotator [--debug|--info] [--cpus=<int>] [--ref-ver=<str>]
        [--normalize-vcf] [--dest-dir=<path>] <data_dir_path> <fa_path>
        <vcf_path>...
    vcline -h|--help
    vcline --version

Commands:
    init                    Create a config YAML template
    run                     Run the analytical pipeline
    download-resources      Download and process resource data
    download-funcotator-data
                            Download Funcotator data sources
    download-gnomad-af-only-vcf
                            Download gnomAD VCF and process it into AF-only VCF
    write-af-only-vcf       Extract and write only AF from VCF INFO
    create-interval-list    Create an interval_list from BED

Options:
    -h, --help              Print help and exit
    --version               Print version and exit
    --debug, --info         Execute a command with debug|info messages
    --yml=<path>            Specify a config YAML path [default: vcline.yml]
    --cpus=<int>            Limit CPU cores used
    --workers=<int>         Specify the maximum number of workers [default: 2]
    --split-intervals       Split evaluation intervals
    --skip-cleaning         Skip incomlete file removal when a task fails
    --ref-dir=<path>        Specify a reference directory path
    --dest-dir=<path>       Specify a destination directory path [default: .]
    --without-gnomad        Skip downloading gnomAD VCF (>200GB)
    --src-url=<url>         Specify a source URL
    --src-path=<path>       Specify a source path
    --ref-ver=<str>         Specify a reference version [default: hg38]
    --normalize-vcf         Normalize VCF files

Args:
    <fa_path>               Path to an reference FASTA file
                            (The index and sequence dictionary are required.)
    <bed_path>              Path to a BED file
    <data_dir_path>         Path to a Funcotator data source directory
    <vcf_path>              Path to a VCF file
"""

import logging
import os
from datetime import datetime
from math import floor
from pathlib import Path

from docopt import docopt
from psutil import cpu_count, virtual_memory

from .. import __version__
from ..task.funcotator import RunFuncotator
from ..task.pipeline import RunAnalyticalPipeline
from ..task.resource import (CreateIntervalListWithBED,
                             DownloadAndConvertVCFsIntoPassingAfOnlyVCF,
                             DownloadFuncotatorDataSources,
                             DownloadResourceFile, WritePassingAfOnlyVCF)
from .util import (build_luigi_tasks, convert_url_to_dest_file_path,
                   fetch_executable, load_default_url_dict, print_log,
                   render_template, write_config_yml)


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
    logger = logging.getLogger(__name__)
    logger.debug(f'args:{os.linesep}{args}')
    if args['init']:
        write_config_yml(path=args['--yml'])
    elif args['run']:
        _run_analytical_pipeline(
            config_yml_path=args['--yml'], ref_dir_path=args['--ref-dir'],
            dest_dir_path=args['--dest-dir'], max_n_cpu=args['--cpus'],
            max_n_worker=args['--workers'],
            split_intervals=args['--split-intervals'],
            skip_cleaning=args['--skip-cleaning'], log_level=log_level
        )
    else:
        dest_dir_path = str(Path(args['--dest-dir']).resolve())
        n_cpu = int(args['--cpus'] or cpu_count())
        if args['download-resources']:
            url_dict = load_default_url_dict()
            cmds = {
                c: fetch_executable(c)
                for c in ['wget', 'pbzip2', 'bgzip', 'gatk']
            }
            build_luigi_tasks(
                tasks=[
                    *[
                        DownloadResourceFile(
                            src_url=v, dest_dir_path=dest_dir_path,
                            n_cpu=n_cpu, wget=cmds['wget'],
                            pbzip2=cmds['pbzip2'], bgzip=cmds['bgzip']
                        ) for k, v in url_dict.items()
                        if k not in {'gnomad_vcf', 'genomesize_xml', 'kmer_fa'}
                    ],
                    *[
                        DownloadResourceFile(
                            src_url=v,
                            dest_file_path=(
                                convert_url_to_dest_file_path(
                                    url=v, dest_dir_path=dest_dir_path
                                ) + '.bz2'
                            ),
                            n_cpu=n_cpu, wget=cmds['wget'],
                            pbzip2=cmds['pbzip2'], bgzip=cmds['bgzip']
                        ) for k, v in url_dict.items()
                        if k in {'genomesize_xml', 'kmer_fa'}
                    ],
                    DownloadFuncotatorDataSources(
                        dest_dir_path=dest_dir_path, n_cpu=n_cpu,
                        gatk=cmds['gatk']
                    ),
                    *(
                        list() if args['--without-gnomad'] else [
                            DownloadAndConvertVCFsIntoPassingAfOnlyVCF(
                                src_url=url_dict['gnomad_vcf'],
                                dest_dir_path=dest_dir_path, n_cpu=n_cpu,
                                wget=cmds['wget'], bgzip=cmds['bgzip']
                            )
                        ]
                    )
                ],
                log_level=log_level
            )
        elif args['download-funcotator-data']:
            build_luigi_tasks(
                tasks=[
                    DownloadFuncotatorDataSources(
                        dest_dir_path=dest_dir_path, n_cpu=n_cpu,
                        gatk=fetch_executable('gatk')
                    )
                ],
                log_level=log_level
            )
        elif args['download-gnomad-af-only-vcf']:
            build_luigi_tasks(
                tasks=[
                    DownloadAndConvertVCFsIntoPassingAfOnlyVCF(
                        src_url=load_default_url_dict()['gnomad_vcf'],
                        dest_dir_path=dest_dir_path, n_cpu=n_cpu,
                        **{c: fetch_executable(c) for c in ['wget', 'bgzip']}
                    )
                ],
                log_level=log_level
            )
        elif args['write-af-only-vcf']:
            build_luigi_tasks(
                tasks=[
                    WritePassingAfOnlyVCF(
                        src_path=(
                            str(Path(args['--src-path']).resolve())
                            if args['--src-path'] else ''
                        ),
                        dest_dir_path=dest_dir_path, n_cpu=n_cpu,
                        **(
                            {'src_url': args['--src-url']}
                            if args['--src-url'] else dict()
                        ),
                        **{c: fetch_executable(c) for c in ['curl', 'bgzip']}
                    )
                ],
                log_level=log_level
            )
        elif args['create-interval-list']:
            build_luigi_tasks(
                tasks=[
                    CreateIntervalListWithBED(
                        fa_path=str(Path(args['<fa_path>']).resolve()),
                        bed_path=str(Path(args['<bed_path>']).resolve()),
                        dest_dir_path=dest_dir_path, n_cpu=n_cpu,
                        gatk=fetch_executable('gatk')
                    )
                ],
                log_level=log_level
            )
        elif args['funcotator']:
            n_vcf = len(args['<vcf_path>'])
            n_worker = min(n_cpu, n_vcf)
            kwargs = {
                'data_src_dir_path':
                str(Path(args['<data_dir_path>']).resolve()),
                'fa_path': str(Path(args['<fa_path>']).resolve()),
                'ref_version': args['--ref-ver'],
                'dest_dir_path': dest_dir_path,
                'normalize_vcf': args['--normalize-vcf'],
                'n_cpu': (floor(n_cpu / n_vcf) if n_cpu > n_vcf else 1),
                'memory_mb':
                int(virtual_memory().total / 1024 / 1024 / n_worker),
                **{c: fetch_executable(c) for c in ['gatk', 'bcftools']}
            }
            build_luigi_tasks(
                tasks=[
                    RunFuncotator(
                        input_vcf_path=str(Path(p).resolve()), **kwargs
                    ) for p in args['<vcf_path>']
                ],
                workers=n_worker, log_level=log_level
            )


def _run_analytical_pipeline(config_yml_path, dest_dir_path='.',
                             ref_dir_path=None, max_n_cpu=None,
                             max_n_worker=None, split_intervals=False,
                             skip_cleaning=False, log_level='WARNING'):
    dest_dir = Path(dest_dir_path).resolve()
    assert dest_dir.is_dir()
    assert Path(config_yml_path).is_file()
    log_dir = dest_dir.joinpath('log')
    if not log_dir.is_dir():
        print_log(f'Make a directory:\t{log_dir}')
        log_dir.mkdir()
    file_log_level = 'DEBUG' if log_level in {'DEBUG', 'INFO'} else 'INFO'
    log_txt_path = str(
        log_dir.joinpath(
            'luigi.{0}.{1}.log.txt'.format(
                datetime.now().strftime('%Y%m%d_%H%M%S'), file_log_level
            )
        )
    )
    log_cfg_path = str(log_dir.joinpath('luigi.log.cfg'))
    render_template(
        template='{}.j2'.format(Path(log_cfg_path).name),
        data={
            'console_log_level': log_level, 'file_log_level': file_log_level,
            'log_txt_path': log_txt_path
        },
        output_path=log_cfg_path
    )
    n_worker = int(max_n_worker or max_n_cpu or cpu_count())
    print_log(f'Run the analytical pipeline:\t{dest_dir}')
    build_luigi_tasks(
        tasks=[
            RunAnalyticalPipeline(
                config_yml_path=config_yml_path, dest_dir_path=str(dest_dir),
                log_dir_path=str(log_dir),
                ref_dir_path=str(
                    Path(ref_dir_path).resolve() if ref_dir_path
                    else dest_dir.joinpath('ref')
                ),
                n_cpu_per_worker=max(
                    1, floor((max_n_cpu or cpu_count()) / n_worker)
                ),
                memory_mb_per_worker=int(
                    virtual_memory().total / 1024 / 1024 / n_worker
                ),
                split_intervals=split_intervals, skip_cleaning=skip_cleaning,
                log_level=log_level
            )
        ],
        workers=n_worker, log_level=log_level, logging_conf_file=log_cfg_path
    )
