#!/usr/bin/env python

import logging
import os
from math import floor
from multiprocessing import cpu_count
from pathlib import Path
from pprint import pformat

import luigi

from ..task.align import AlignReads
from .util import (fetch_executable, is_url, print_log, read_yml,
                   render_template)


def run_analytical_pipeline(config_yml_path, work_dir_path=None,
                            ref_dir_path=None, max_n_cpu=None,
                            max_n_worker=None, log_level='WARNING'):
    logger = logging.getLogger(__name__)
    work_dir = Path(work_dir_path or '.').resolve()
    print_log(f'Run the analytical pipeline:\t{work_dir}')
    assert work_dir.is_dir()
    cf = _read_config_yml(config_yml_path=config_yml_path)
    logger.debug('cf:' + os.linesep + pformat(cf))

    log_dir = work_dir.joinpath('log')
    dirs = {
        **{
            f'{k}_dir_path': str(work_dir.joinpath(k))
            for k in ['trim', 'align']
        },
        'ref_dir_path':
        str(Path(ref_dir_path or str(work_dir.joinpath('ref'))).resolve()),
        'log_dir_path': str(log_dir)
    }
    n_worker = max(1, int(max_n_worker or 1))
    n_cpu_per_worker = max(1, floor(int(max_n_cpu or cpu_count()) / n_worker))
    ref_fa_list = [{'src': u, 'is_url': is_url(u)} for u in cf['ref_fa']]
    common_params = {
        **{
            c: fetch_executable(c) for c in [
                'samtools', 'bwa', 'fastqc', 'cutadapt', 'trim_galore', 'pigz',
                'pbzip2', 'cat', 'curl'
            ]
        },
        **dirs,
        'n_cpu_per_worker': n_cpu_per_worker
    }
    logger.debug('common_params:' + os.linesep + pformat(common_params))

    for p in dirs.values():
        d = Path(p)
        if not d.is_dir():
            print_log(f'Make a directory:\t{p}')
            d.mkdir()
    luigi_log_cfg_path = str(log_dir.joinpath('luigi.log.cfg'))
    render_template(
        template='{}.j2'.format(Path(luigi_log_cfg_path).name),
        data={
            'level': log_level,
            'filename': str(log_dir.joinpath('luigi.log.txt'))
        },
        output_path=luigi_log_cfg_path
    )

    for r in cf['runs']:
        luigi.build(
            [
                AlignReads(
                    fq_paths=[str(Path(p).resolve()) for p in r[k]['fq']],
                    ref_fa_list=ref_fa_list, p=common_params
                ) for k in ['foreground', 'background'] if r[k].get('fq')
            ],
            workers=n_worker, local_scheduler=True, log_level=log_level,
            logging_conf_file=luigi_log_cfg_path
        )


def _read_config_yml(config_yml_path):
    cf = read_yml(path=str(Path(config_yml_path).resolve()))
    assert isinstance(cf, dict)
    for k in ['ref_fa', 'runs']:
        assert cf.get(k)
        assert isinstance(cf[k], list)
    assert _has_unique_elements(cf['ref_fa'])
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
                assert _has_unique_elements(r[k]['fq'])
                assert len(r[k]['fq']) <= 2
            else:
                assert isinstance(r[k]['sam'], str)
    assert _has_unique_elements([r['id'] for r in cf['runs']])
    return cf


def _has_unique_elements(elements):
    return len(set(elements)) == len(tuple(elements))
