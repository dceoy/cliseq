#!/usr/bin/env python
"""
Variant Calling Pipeline for Clinical Sequencing

Usage:
    vcline init [--debug|--info] [--yml=<path>]
    vcline download [--debug|--info] [--cpus=<int>] [--workers=<int>]
        [--skip-cleaning] [--print-subprocesses] [--use-bwa-mem2]
        [--snpeff|--funcotator|--vep] [--dest-dir=<path>]
    vcline run [--debug|--info] [--yml=<path>] [--cpus=<int>]
        [--workers=<int>] [--skip-cleaning] [--print-subprocesses]
        [--use-bwa-mem2] [--snpeff|--funcotator|--vep] [--dest-dir=<path>]
    vcline write-af-only-vcf [--debug|--info] [--cpus=<int>]
        [--src-path=<path>|--src-url=<url>] [--dest-dir=<path>]
    vcline -h|--help
    vcline --version

Commands:
    init                    Create a config YAML template
    download                Download and preprocess hg38 resources
    run                     Run the analytical pipeline
    write-af-only-vcf       Extract and write only AF from VCF INFO

Options:
    -h, --help              Print help and exit
    --version               Print version and exit
    --debug, --info         Execute a command with debug|info messages
    --yml=<path>            Specify a config YAML path [default: vcline.yml]
    --cpus=<int>            Limit CPU cores used
    --workers=<int>         Specify the maximum number of workers [default: 1]
    --skip-cleaning         Skip incomlete file removal when a task fails
    --print-subprocesses    Print STDOUT/STDERR outputs from subprocesses
    --use-bwa-mem2          Use BWA-MEM2 for read alignment
    --snpeff, --funotator, --vep
                            Select only one of SnpEff, Funcotator, and VEP
    --dest-dir=<path>       Specify a destination directory path [default: .]
    --src-path=<path>       Specify a source path
    --src-url=<url>         Specify a source URL
"""

import logging
import os
from pathlib import Path

from docopt import docopt
from ftarc.cli.util import build_luigi_tasks, fetch_executable, print_log
from psutil import cpu_count, virtual_memory
from vanqc.task.gatk import DownloadFuncotatorDataSources
from vanqc.task.snpeff import DownloadSnpeffDataSources
from vanqc.task.vep import DownloadEnsemblVepCache

from .. import __version__
from ..task.downloader import PreprocessResources, WritePassingAfOnlyVcf
from .pipeline import run_analytical_pipeline
from .util import load_default_dict, write_config_yml


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
    print_log(f'Start the workflow of vcline {__version__}')
    if args['init']:
        write_config_yml(path=args['--yml'])
    elif args['run']:
        run_analytical_pipeline(
            config_yml_path=args['--yml'], dest_dir_path=args['--dest-dir'],
            max_n_cpu=args['--cpus'], max_n_worker=args['--workers'],
            skip_cleaning=args['--skip-cleaning'],
            print_subprocesses=args['--print-subprocesses'],
            console_log_level=log_level, use_bwa_mem2=args['--use-bwa-mem2']
        )
    else:
        n_cpu = int(args['--cpus'] or cpu_count())
        memory_mb = virtual_memory().total / 1024 / 1024 / 2
        sh_config = {
            'log_dir_path': None,
            'remove_if_failed': (not args['--skip-cleaning']), 'quiet': False,
            'executable': fetch_executable('bash')
        }
        if args['download']:
            url_dict = load_default_dict(stem='urls')
            command_dict = {
                'bwa': fetch_executable(
                    'bwa-mem2' if args['--use-bwa-mem2'] else 'bwa'
                ),
                **{
                    c: fetch_executable(c) for c in [
                        'wget', 'pbzip2', 'bgzip', 'pigz', 'samtools', 'tabix',
                        'gatk'
                    ]
                }
            }
            anns = (
                {k for k in ['snpeff', 'funcotator', 'vep'] if args[f'--{k}']}
                or {'snpeff', 'funcotator', 'vep'}
            )
            common_kwargs = {
                'dest_dir_path': args['--dest-dir'], 'sh_config': sh_config
            }
            build_luigi_tasks(
                tasks=[
                    PreprocessResources(
                        src_url_dict=url_dict, **command_dict, n_cpu=n_cpu,
                        memory_mb=memory_mb,
                        use_bwa_mem2=args['--use-bwa-mem2'], **common_kwargs
                    ),
                    *(
                        [
                            DownloadSnpeffDataSources(
                                snpeff=fetch_executable('snpEff'),
                                genome_version='GRCh38', memory_mb=memory_mb,
                                **common_kwargs
                            )
                        ] if 'snpeff' in anns else list()
                    ),
                    *(
                        [
                            DownloadFuncotatorDataSources(
                                gatk=command_dict['gatk'], n_cpu=n_cpu,
                                memory_mb=memory_mb, **common_kwargs
                            )
                        ] if 'funcotator' in anns else list()
                    ),
                    *(
                        [
                            DownloadEnsemblVepCache(
                                genome_version='GRCh38',
                                vep=fetch_executable('vep'),
                                wget=command_dict['wget'], **common_kwargs
                            )
                        ] if 'vep' in anns else list()
                    )
                ],
                log_level=log_level
            )
        elif args['write-af-only-vcf']:
            dest_dir_path = str(Path(args['--dest-dir']).resolve())
            build_luigi_tasks(
                tasks=[
                    WritePassingAfOnlyVcf(
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
