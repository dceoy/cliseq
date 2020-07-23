#!/usr/bin/env python

import logging
import os
import re
from datetime import datetime
from itertools import chain, product
from math import floor
from pathlib import Path
from pprint import pformat

import luigi
import yaml
from psutil import cpu_count, virtual_memory

from ..cli.util import (fetch_executable, load_default_dict, parse_cram_id,
                        parse_fq_id, print_log, read_yml, render_template)
from ..task.pipeline import (PrepareCRAMsMatched, PreprocessResources,
                             PrintEnvVersions, RunVariantCaller)


def build_luigi_tasks(*args, **kwargs):
    r = luigi.build(
        *args,
        **{
            k: v for k, v in kwargs.items() if (
                k not in {'logging_conf_file', 'hide_summary'}
                or (k == 'logging_conf_file' and v)
            )
        },
        local_scheduler=True, detailed_summary=True
    )
    if not kwargs.get('hide_summary'):
        print(
            os.linesep
            + os.linesep.join(['Execution summary:', r.summary_text, str(r)])
        )


def run_analytical_pipeline(config_yml_path, dest_dir_path=None,
                            ref_dir_path=None, max_n_cpu=None,
                            max_n_worker=None, skip_cleaning=False,
                            print_subprocesses=False,
                            console_log_level='WARNING',
                            file_log_level='DEBUG', only_preprocessing=False,
                            use_bwa_mem2=True):
    logger = logging.getLogger(__name__)
    logger.info(f'config_yml_path:\t{config_yml_path}')
    config = _read_config_yml(
        path=config_yml_path, only_preprocessing=only_preprocessing
    )
    runs = config.get('runs')
    logger.info(f'dest_dir_path:\t{dest_dir_path}')
    dest_dir = Path(dest_dir_path).resolve()
    log_dir = dest_dir.joinpath('log')

    read_alignment = bool(
        runs and [r for r in runs if [v for v in r.values() if v.get('fq')]]
    )
    logger.debug(f'read_alignment:\t{read_alignment}')
    adapter_removal = (
        read_alignment and (
            config['adapter_removal'] if 'adapter_removal' in config
            else True
        )
    )
    logger.debug(f'adapter_removal:\t{adapter_removal}')

    default_dict = load_default_dict(stem='vcline')
    callers = (
        list(
            chain.from_iterable([
                [
                    f'{k}.{c}' for c, b in v.items()
                    if b and c in default_dict['callers'][k]
                ] for k, v in config['callers'].items()
                if default_dict['callers'].get(k)
            ])
        ) if 'callers' in config else list()
    )
    logger.debug('callers:' + os.linesep + pformat(callers))

    annotators = (
        [k for k in default_dict['annotators'] if config['annotators'].get(k)]
        if 'annotators' in config else list()
    )
    logger.debug('annotators:' + os.linesep + pformat(annotators))

    command_dict = {
        c: fetch_executable(c) for c in {
            'bgzip', 'gatk', 'java', 'pbzip2', 'pigz', 'samtools', 'tabix',
            *({'bcftools'} if only_preprocessing or callers else set()),
            *(
                {'cutadapt', 'fastqc', 'trim_galore'}
                if adapter_removal else set()
            ),
            *(
                {'python2', 'configureStrelkaSomaticWorkflow.py'}
                if 'somatic_snv_indel.strelka' in callers else set()
            ),
            *({'python3'} if 'germline_snv_indel.gatk' in callers else set()),
            *(
                {'python2', 'configureStrelkaGermlineWorkflow.py'}
                if 'germline_snv_indel.strelka' in callers else set()
            ),
            *(
                {'python2', 'configManta.py'}
                if 'somatic_sv.manta' in callers else set()
            ),
            *(
                {'bedtools'}
                if only_preprocessing or 'somatic_sv.delly' in callers
                else set()
            ),
            *({'delly'} if 'somatic_sv.delly' in callers else set()),
            *({'R'} if 'somatic_cnv.gatk' in callers else set()),
            *({'snpEff'} if 'snpeff' in annotators else set()),
            *(
                {'msisensor'}
                if only_preprocessing or 'somatic_msi.msisensor' in callers
                else set()
            )
        }
    }
    if only_preprocessing or read_alignment:
        command_dict['bwa'] = fetch_executable(
            'bwa-mem2' if use_bwa_mem2 else 'bwa'
        )
    logger.debug('command_dict:' + os.linesep + pformat(command_dict))

    n_cpu = cpu_count()
    n_worker = min(
        int(max_n_worker or max_n_cpu or n_cpu),
        (((len(callers) if callers else 2) * len(runs)) if runs else 16)
    )
    n_cpu_per_worker = max(1, floor((max_n_cpu or n_cpu) / n_worker))
    memory_mb = virtual_memory().total / 1024 / 1024 / 2
    memory_mb_per_worker = int(memory_mb / n_worker)
    cf_dict = {
        'log_dir_path': str(log_dir),
        'ref_dir_path':
        (str(Path(ref_dir_path).resolve()) if ref_dir_path else None),
        'n_worker': n_worker, 'memory_mb_per_worker': memory_mb_per_worker,
        'n_cpu_per_worker': n_cpu_per_worker,
        'gatk_java_options': ' '.join([
            '-Dsamjdk.compression_level=5',
            '-Dsamjdk.use_async_io_read_samtools=true',
            '-Dsamjdk.use_async_io_write_samtools=true',
            '-Dsamjdk.use_async_io_write_tribble=false',
            f'-Xmx{memory_mb_per_worker:d}m',
            '-XX:+UseParallelGC',
            f'-XX:ParallelGCThreads={n_cpu_per_worker}'
        ]),
        'ref_version': (config.get('reference_version') or 'hg38'),
        'exome': bool(config.get('exome')),
        'adapter_removal': adapter_removal,
        'save_memory': (memory_mb_per_worker < 8 * 1024),
        'remove_if_failed': (not skip_cleaning),
        'quiet': (not print_subprocesses),
        **{
            (k.replace('/', '_') + '_dir_path'): str(dest_dir.joinpath(k))
            for k in {
                'trim', 'align', 'qc', 'postproc',
                *chain.from_iterable([
                    [f'{k}/{c}' for c in v.keys()]
                    for k, v in default_dict['callers'].items()
                ])
            }
        },
        **command_dict
    }
    logger.debug('cf_dict:' + os.linesep + pformat(cf_dict))

    if only_preprocessing:
        resource_keys = {
            'ref_fa', 'dbsnp_vcf', 'mills_indel_vcf', 'known_indel_vcf',
            'evaluation_interval', 'hapmap_vcf', 'gnomad_vcf', 'cnv_blacklist'
        }
    elif callers:
        resource_keys = {
            'ref_fa', 'dbsnp_vcf', 'mills_indel_vcf', 'known_indel_vcf',
            'evaluation_interval', 'hapmap_vcf', 'gnomad_vcf', 'cnv_blacklist',
            'funcotator_somatic_tar', 'funcotator_germline_tar',
            'snpeff_config'
        }
    else:
        resource_keys = {
            'ref_fa', 'dbsnp_vcf', 'mills_indel_vcf', 'known_indel_vcf'
        }
    resource_path_dict = _resolve_input_file_paths(
        path_dict={
            k: v for k, v in config['resources'].items() if k in resource_keys
        }
    )
    logger.debug(
        'resource_path_dict:' + os.linesep + pformat(resource_path_dict)
    )

    sample_dict_list = (
        [
            {**_determine_input_samples(run_dict=r), 'priority': p} for p, r
            in zip([i * 1000 for i in range(1, (len(runs) + 1))[::-1]], runs)
        ] if runs else list()
    )
    logger.debug('sample_dict_list:' + os.linesep + pformat(sample_dict_list))

    if only_preprocessing:
        print_log(f'Run the resources preprocessing:\t{dest_dir}')
        print(
            yaml.dump([
                {'workers': n_worker}, {'resources': resource_path_dict}
            ])
        )
    else:
        print_log(
            (
                'Run the analytical pipeline' if callers
                else 'Prepare CRAM files'
            ) + f'\t{dest_dir}'
        )
        print(
            yaml.dump([
                {'workers': n_worker}, {'runs': len(runs)},
                {'adapter_removal': adapter_removal},
                {'read_alignment': read_alignment},
                {'callers': callers}, {'annotators': annotators},
                {
                    'samples': [
                        dict(zip(['tumor', 'normal'], d['sample_names']))
                        for d in sample_dict_list
                    ]
                }
            ])
        )

    for d in [dest_dir, log_dir]:
        if not d.is_dir():
            print_log(f'Make a directory:\t{d}')
            d.mkdir()
    log_cfg_path = str(log_dir.joinpath('luigi.log.cfg'))
    render_template(
        template=(Path(log_cfg_path).name + '.j2'),
        data={
            'console_log_level': console_log_level,
            'file_log_level': file_log_level,
            'log_txt_path': str(
                log_dir.joinpath(
                    'luigi.{0}.{1}.log.txt'.format(
                        file_log_level,
                        datetime.now().strftime('%Y%m%d_%H%M%S')
                    )
                )
            )
        },
        output_path=log_cfg_path
    )

    build_luigi_tasks(
        tasks=[
            PrintEnvVersions(
                log_dir_path=str(log_dir),
                command_paths=list(command_dict.values())
            )
        ],
        workers=1, log_level=console_log_level, logging_conf_file=log_cfg_path,
        hide_summary=True
    )
    if only_preprocessing:
        build_luigi_tasks(
            tasks=[PreprocessResources(**resource_path_dict, cf=cf_dict)],
            workers=n_worker, log_level=console_log_level,
            logging_conf_file=log_cfg_path
        )
    elif callers:
        build_luigi_tasks(
            tasks=[
                RunVariantCaller(
                    **d, **resource_path_dict, cf=cf_dict, caller=c,
                    annotators=annotators
                ) for d, c in product(sample_dict_list, callers)
            ],
            workers=n_worker, log_level=console_log_level,
            logging_conf_file=log_cfg_path
        )
    else:
        build_luigi_tasks(
            tasks=[
                PrepareCRAMsMatched(**d, **resource_path_dict, cf=cf_dict)
                for d in sample_dict_list
            ],
            workers=n_worker, log_level=console_log_level,
            logging_conf_file=log_cfg_path
        )


def _read_config_yml(path, only_preprocessing=False):
    config = read_yml(path=Path(path).resolve())
    assert (isinstance(config, dict) and config.get('resources')), config
    assert isinstance(config['resources'], dict), config['resources']
    for k in ['ref_fa', 'dbsnp_vcf', 'mills_indel_vcf', 'known_indel_vcf',
              'hapmap_vcf', 'gnomad_vcf', 'evaluation_interval',
              'funcotator_germline_tar', 'funcotator_somatic_tar',
              'snpeff_config']:
        v = config['resources'].get(k)
        if k == 'ref_fa' and isinstance(v, list) and v:
            assert _has_unique_elements(v), k
            for s in v:
                assert isinstance(s, str), k
        elif v:
            assert isinstance(v, str), k
    if not only_preprocessing:
        assert config.get('runs'), config
        assert isinstance(config['runs'], list), config['runs']
        for r in config['runs']:
            assert isinstance(r, dict), r
            assert set(r.keys()).intersection({'tumor', 'normal'}), r
            for t in ['tumor', 'normal']:
                s = r[t]
                assert isinstance(s, dict), s
                assert (s.get('fq') or s.get('cram')), s
                if s.get('fq'):
                    assert isinstance(s['fq'], list), s
                    assert _has_unique_elements(s['fq']), s
                    assert (len(s['fq']) <= 2), s
                    for p in s['fq']:
                        assert p.endswith(('.gz', '.bz2')), p
                else:
                    assert isinstance(s['cram'], str), s
                    assert s['cram'].endswith('.cram'), s
                if s.get('read_group'):
                    assert isinstance(s['read_group'], dict), s
                    for k, v in s['read_group'].items():
                        assert re.fullmatch(r'[A-Z]{2}', k), k
                        assert isinstance(v, str), k
                if s.get('sample_name'):
                    assert isinstance(s['sample_name'], str), s
    return config


def _has_unique_elements(elements):
    return len(set(elements)) == len(tuple(elements))


def _resolve_file_path(path):
    p = Path(path).resolve()
    assert p.is_file(), f'file not found: {p}'
    return str(p)


def _resolve_input_file_paths(path_list=None, path_dict=None):
    if path_list:
        return [_resolve_file_path(s) for s in path_list]
    elif path_dict:
        new_dict = dict()
        for k, v in path_dict.items():
            if isinstance(v, str):
                new_dict[f'{k}_path'] = _resolve_file_path(v)
            elif v:
                new_dict[f'{k}_paths'] = [
                    _resolve_file_path(s) for s in v
                ]
        return new_dict


def _determine_input_samples(run_dict):
    tn = [run_dict[i] for i in ['tumor', 'normal']]
    cram_paths = [d.get('cram') for d in tn]
    if all(cram_paths):
        cram_list = _resolve_input_file_paths(path_list=cram_paths)
        return {
            'fq_list': list(), 'read_groups': list(),
            'cram_list': cram_list,
            'sample_names': [
                (d.get('sample_name') or parse_cram_id(cram_path=d['cram']))
                for d in tn
            ]
        }
    else:
        return {
            'fq_list':
            [list(_resolve_input_file_paths(path_list=d['fq'])) for d in tn],
            'read_groups': [(d.get('read_group') or dict()) for d in tn],
            'cram_list': list(),
            'sample_names': [
                (
                    (d.get('read_group') or dict()).get('SM')
                    or parse_fq_id(fq_path=d['fq'][0])
                ) for d in tn
            ]
        }
