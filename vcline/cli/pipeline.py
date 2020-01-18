#!/usr/bin/env python

import logging
import os
from multiprocessing import cpu_count
from pathlib import Path
from pprint import pformat

import luigi

from ..task.ref import PrepareReferences
from .util import fetch_executable, print_log, read_yml


def run_analytical_pipeline(config_yml_path, work_dir_path=None,
                            ref_dir_path=None, n_cpu=None):
    logger = logging.getLogger(__name__)
    work_dir = Path(work_dir_path or '.').resolve()
    print_log('Run the analytical pipeline:\t{}'.format(work_dir))
    assert work_dir.is_dir()

    config_yml_abspath = str(Path(config_yml_path).resolve())
    cf = read_yml(path=config_yml_abspath)
    logger.debug('cf:' + os.linesep + pformat(cf))
    assert 'ref_fa' in cf
    assert set(cf['ref_fa']).intersection({'urls', 'paths'})

    param_dict = {
        **cf,
        **{
            c: fetch_executable(c) for c in [
                'samtools', 'bwa', 'fastqc', 'cutadapt', 'trim_galore', 'pigz'
            ]
        },
        'ref_dir_path': str(
            Path(ref_dir_path).resolve()
            if ref_dir_path else work_dir.joinpath('ref')
        ),
        'n_cpu': int(n_cpu or cpu_count())
    }
    logger.debug('param_dict:' + os.linesep + pformat(param_dict))

    luigi.build(
        [
            PrepareReferences(**{
                k: v for k, v in param_dict.items() if k in [
                    'ref_fa', 'ref_dir_path', 'samtools', 'bwa'
                ]
            })
        ],
        workers=2, local_scheduler=True
    )
