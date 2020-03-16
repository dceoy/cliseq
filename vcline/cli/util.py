#!/usr/bin/env python

import logging
import os
import re
from pathlib import Path
from pprint import pformat

import luigi
import yaml
from jinja2 import Environment, FileSystemLoader


def print_log(message):
    logger = logging.getLogger(__name__)
    logger.debug(message)
    print(f'>>\t{message}', flush=True)


def build_luigi_tasks(*args, **kwargs):
    print_log(
        'Task execution summary:' + os.linesep + luigi.build(
            *args, **kwargs, local_scheduler=True, detailed_summary=True
        ).summary_text
    )


def parse_fq_id(fq_path):
    return (
        re.sub(
            r'([\._]read[12][\._]|[\._]r[12][\._]|\.fq\.|\.fastq\.).+$', '',
            Path(fq_path).name, flags=re.IGNORECASE
        ) or Path(Path(fq_path).stem).stem
    )


def create_matched_id(tumor_name, normal_name):
    fragments = [Path(n).stem.split('.') for n in [tumor_name, normal_name]]
    if fragments[0][-1] != fragments[1][-1]:
        somatic_id = '.'.join(['.'.join(f) for f in fragments])
    else:
        n_common = 0
        for i in range(1, min([len(f) for f in fragments])):
            if fragments[0][-i] == fragments[1][-i]:
                n_common += 1
            else:
                break
        somatic_id = '.'.join(fragments[0][:-n_common] + fragments[1])
    return somatic_id


def fetch_executable(cmd):
    executables = [
        cp for cp in [
            str(Path(p).joinpath(cmd))
            for p in os.environ['PATH'].split(os.pathsep)
        ] if os.access(cp, os.X_OK)
    ]
    if executables:
        return executables[0]
    else:
        raise RuntimeError(f'command not found: {cmd}')


def read_config_yml(config_yml_path):
    config = _read_yml(path=str(Path(config_yml_path).resolve()))
    assert isinstance(config, dict)
    for k in ['references', 'runs']:
        assert config.get(k)
    assert isinstance(config['references'], dict)
    assert {
        'ref_fa', 'dbsnp_vcf', 'mills_indel_vcf', 'known_indel_vcf',
        'hapmap_vcf', 'gnomad_vcf', 'evaluation_interval',
        'funcotator_germline_tar', 'funcotator_somatic_tar'
    }.issubset(set(config['references'].keys()))
    for k in ['ref_fa', 'dbsnp_vcf', 'mills_indel_vcf', 'known_indel_vcf',
              'hapmap_vcf', 'gnomad_vcf', 'evaluation_interval',
              'funcotator_germline_tar', 'funcotator_somatic_tar']:
        v = config['references'].get(k)
        if k == 'ref_fa' and isinstance(v, list):
            assert _has_unique_elements(v)
            for s in v:
                assert isinstance(s, str)
        else:
            assert isinstance(v, str)
    assert isinstance(config['runs'], list)
    for r in config['runs']:
        assert isinstance(r, dict)
        assert set(r.keys()).intersection({'tumor', 'normal'})
        for t in ['tumor', 'normal']:
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


def _read_yml(path):
    logger = logging.getLogger(__name__)
    with open(path, 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    logger.debug('YAML data:' + os.linesep + pformat(d))
    return d


def _has_unique_elements(elements):
    return len(set(elements)) == len(tuple(elements))


def render_template(template, data, output_path):
    print_log(
        ('Overwrite' if Path(output_path).exists() else 'Render')
        + f' a file:\t{output_path}'
    )
    with open(output_path, 'w') as f:
        f.write(
            Environment(
                loader=FileSystemLoader(
                    str(Path(__file__).parent.joinpath('../template')),
                    encoding='utf8'
                )
            ).get_template(template).render(data) + os.linesep
        )


def load_default_url_dict(include_gnomad=False):
    d = _read_yml(
        path=str(Path(__file__).parent.parent.joinpath('static/urls.yml'))
    )
    return (
        d if include_gnomad
        else {k: v for k, v in d.items() if k != 'gnomad_vcf'}
    )


def generate_gnomad_chrom_vcf_urls():
    url_suffix = re.sub(
        r'\.vcf\.bgz$', '',
        load_default_url_dict(include_gnomad=True)['gnomad_vcf']
    )
    for c in [*range(1, 23), 'X', 'Y']:
        yield f'{url_suffix}.chr{c}.vcf.bgz'
