#!/usr/bin/env python

import logging
import os
from multiprocessing import cpu_count
from pathlib import Path
from pprint import pformat

import luigi

from ..task.ref import CreateBWAIndexes, CreateFASTAIndex
from ..task.trim import TrimAdapters
from .util import fetch_executable, print_log, read_yml


def run_analytical_pipeline(config_yml_path, work_dir_path=None,
                            ref_dir_path=None, n_cpu=None,
                            log_level='WARNING'):
    logger = logging.getLogger(__name__)
    work_dir = Path(work_dir_path or '.').resolve()
    print_log('Run the analytical pipeline:\t{}'.format(work_dir))
    assert work_dir.is_dir()
    cf = _read_config_yml(config_yml_path=config_yml_path)
    logger.debug('cf:' + os.linesep + pformat(cf))
    common_params = {
        'ref_fa': cf['ref_fa'],
        **{
            c: fetch_executable(c) for c in [
                'samtools', 'bwa', 'fastqc', 'cutadapt', 'trim_galore',
                'pigz', 'pbzip2'
            ]
        },
        **{
            f'{k}_dir_path': str(work_dir.joinpath(k))
            for k in ['log', 'trim', 'align', 'call']
        },
        'ref_dir_path': str(
            Path(ref_dir_path).resolve()
            if ref_dir_path else work_dir.joinpath('ref')
        ),
        'n_cpu': int(n_cpu or cpu_count())
    }
    fb = ['foreground', 'background']
    for r in cf['runs']:
        params = {
            'run_id': r['id'], **common_params,
            **{f'{k}_fq': (r[k].get('fq') or list()) for k in fb},
            **{f'{k}_sam': r[k].get('sam') for k in fb}
        }
        logger.debug('params:' + os.linesep + pformat(params))
        luigi.build(
            [
                CreateBWAIndexes(**params),
                CreateFASTAIndex(**params),
                *(
                    [TrimAdapters(**params)]
                    if any([params[f'{k}_fq'] for k in fb]) else list()
                )
            ],
            workers=2, local_scheduler=True, log_level=log_level
        )


def _read_config_yml(config_yml_path):
    cf = read_yml(path=str(Path(config_yml_path).resolve()))
    assert isinstance(cf, dict)
    for k in ['ref_fa', 'runs']:
        assert cf.get(k)
        assert isinstance(cf[k], list)
    assert len(cf['ref_fa']) == len(set(cf['ref_fa']))
    for s in cf['ref_fa']:
        assert isinstance(s, str)
    for r in cf['runs']:
        assert isinstance(r, dict)
        assert set(r.keys()).intersection({'id', 'foreground', 'background'})
        for k in ['foreground', 'background']:
            assert isinstance(r[k], dict)
            assert r[k].get('fq') or r[k].get('sam')
            if 'fq' in r[k]:
                assert isinstance(r[k]['fq'], list)
                assert len(r[k]['fq']) == len(set(r[k]['fq']))
                assert len(r[k]['fq']) <= 2
            else:
                assert isinstance(r[k]['sam'], str)
    return cf
