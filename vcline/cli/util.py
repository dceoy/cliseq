#!/usr/bin/env python

import logging
import os
import re
import shutil
from pathlib import Path
from pprint import pformat
from urllib.parse import urlsplit

import luigi
import yaml
from jinja2 import Environment, FileSystemLoader


def write_config_yml(path):
    if Path(path).is_file():
        print_log(f'The file exists:\t{path}')
    else:
        print_log(f'Create a config YAML:\t{path}')
        shutil.copyfile(
            str(Path(__file__).parent.joinpath('../static/vcline.yml')),
            Path(path).resolve()
        )


def convert_url_to_dest_file_path(url, dest_dir_path='.'):
    res = urlsplit(url)
    return str(
        Path(dest_dir_path).joinpath(
            '.'.join((res.netloc.split('.')[0] + res.path).split('/'))
        )
    )


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


def parse_cram_id(cram_path):
    prefix = Path(cram_path).stem
    if '.trim.' not in prefix:
        return prefix
    else:
        t = Path(cram_path).stem.split('.')
        return '.'.join(t[:-(t[::-1].index('trim') + 1)])


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


def read_yml(path):
    logger = logging.getLogger(__name__)
    with open(path, 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    logger.debug('YAML data:' + os.linesep + pformat(d))
    return d


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


def load_default_url_dict():
    return read_yml(
        path=str(Path(__file__).parent.parent.joinpath('static/urls.yml'))
    )
