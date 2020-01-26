#!/usr/bin/env python

import logging
import os
import re
from math import floor
from pathlib import Path
from pprint import pformat

import luigi
from psutil import cpu_count, virtual_memory

from ..task.call import CallVariants
from .util import fetch_executable, print_log, read_yml, render_template


def run_analytical_pipeline(config_yml_path, work_dir_path=None,
                            ref_dir_path=None, max_n_cpu=None,
                            max_n_worker=None, log_level='WARNING'):
    logger = logging.getLogger(__name__)
    work_dir = Path(work_dir_path or '.').resolve()
    print_log(f'Run the analytical pipeline:\t{work_dir}')
    assert work_dir.is_dir()
    config = _read_config_yml(config_yml_path=config_yml_path)
    logger.debug('config:' + os.linesep + pformat(config))

    n_worker = max(1, int(max_n_worker or 1))
    log_dir = work_dir.joinpath('log')
    dirs = {
        **{
            f'{k}_dir_path': str(work_dir.joinpath(k))
            for k in ['trim', 'align', 'call']
        },
        'ref_dir_path':
        str(Path(ref_dir_path or str(work_dir.joinpath('ref'))).resolve()),
        'log_dir_path': str(log_dir)
    }

    total_n_cpu = cpu_count()
    total_memory_mb = virtual_memory().total / 1024 / 1024
    common_config = {
        'memory_mb_per_worker': int(total_memory_mb / n_worker),
        'n_cpu_per_worker':
        max(1, floor(int(max_n_cpu or total_n_cpu) / n_worker)),
        'gatk_java_options':
        '-Dsamjdk.compression_level={0:d} -Xms{1:d}m'.format(
            5, int(total_memory_mb / n_worker / 2)
        ),
        'samtools_memory_per_thread': '{:d}M'.format(
            int(total_memory_mb / total_n_cpu / 20)
        ),
        **{
            c: fetch_executable(c) for c in [
                'bgzip', 'bcftools', 'bwa', 'cat', 'cutadapt', 'fastqc',
                'gatk', 'pbzip2', 'pigz', 'samtools', 'tabix', 'trim_galore'
            ]
        },
        **dirs
    }
    logger.debug('common_config:' + os.linesep + pformat(common_config))
    ref_dict = {
        f'{k}_paths': list(_resolve_input_file_paths(v))
        for k, v in config['references'].items()
    }
    logger.debug('ref_dict:' + os.linesep + pformat(ref_dict))
    log_txt_path = str(log_dir.joinpath('luigi.log.txt'))
    log_cfg_path = str(log_dir.joinpath('luigi.log.cfg'))

    for p in dirs.values():
        d = Path(p)
        if not d.is_dir():
            print_log(f'Make a directory:\t{p}')
            d.mkdir()
    render_template(
        template='{}.j2'.format(Path(log_cfg_path).name),
        data={'log_level': log_level, 'log_txt_path': log_txt_path},
        output_path=log_cfg_path
    )
    luigi.build(
        [
            CallVariants(
                fq_list=[
                    list(_resolve_input_file_paths(paths=r[k]['fq']))
                    for k in ['foreground', 'background'] if r[k].get('fq')
                ],
                read_groups=[
                    (r[k].get('read_group') or dict())
                    for k in ['foreground', 'background']
                ],
                cf=common_config, **ref_dict
            ) for r in config['runs']
        ],
        workers=n_worker, local_scheduler=True, log_level=log_level,
        logging_conf_file=log_cfg_path
    )


def _resolve_input_file_paths(paths):
    for s in paths:
        p = Path(s).resolve()
        assert p.is_file(), f'file not found: {p}'
        yield str(p)


def _read_config_yml(config_yml_path):
    config = read_yml(path=str(Path(config_yml_path).resolve()))
    assert isinstance(config, dict)
    for k in ['references', 'runs']:
        assert config.get(k)
    assert isinstance(config['references'], dict)
    for k in ['ref_fa', 'known_site_vcf']:
        v = config['references'].get(k)
        assert isinstance(v, list)
        assert _has_unique_elements(v)
        for s in v:
            assert isinstance(s, str)
    assert isinstance(config['runs'], list)
    for r in config['runs']:
        assert isinstance(r, dict)
        assert set(r.keys()).intersection({'foreground', 'background'})
        for t in ['foreground', 'background']:
            assert isinstance(r[t], dict)
            assert r[t].get('fq') or r[t].get('sam')
            if r[t].get('fq'):
                assert isinstance(r[t]['fq'], list)
                assert _has_unique_elements(r[t]['fq'])
                assert len(r[t]['fq']) <= 2
            else:
                assert isinstance(r[t]['sam'], str)
            if r[t].get('read_group'):
                assert isinstance(r[k]['read_group'], dict)
                for k, v in r[t]['read_group'].items():
                    assert re.fullmatch(r'[A-Z]{2}', k)
                    assert type(v) not in [list, dict]
    return config


def _has_unique_elements(elements):
    return len(set(elements)) == len(tuple(elements))
