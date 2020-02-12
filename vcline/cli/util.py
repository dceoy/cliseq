#!/usr/bin/env python

import logging
import os
import re
import subprocess
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader


def print_log(message):
    logger = logging.getLogger(__name__)
    logger.debug(message)
    print(f'>>\t{message}', flush=True)


def parse_fq_id(fq_path):
    return (
        re.sub(
            r'([\._]R[12][\._]|\.fq\.|\.fastq\.).+$', '',
            Path(fq_path).name
        ) or Path(Path(fq_path).stem).stem
    )


def create_matched_id(foreground_name, background_name):
    fragments = [
        Path(n).stem.split('.') for n in [foreground_name, background_name]
    ]
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
    with open(path, 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
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
        path=str(Path(__file__).parent.joinpath('../static/urls.yml'))
    )


def run_and_parse_subprocess(args, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, **kwargs):
    logger = logging.getLogger(__name__)
    logger.debug(f'args: {args}')
    with subprocess.Popen(args=args, stdout=stdout, stderr=stderr,
                          **kwargs) as p:
        for line in p.stdout:
            yield line.decode('utf-8')
        outs, errs = p.communicate()
        if p.returncode != 0:
            logger.error(
                'STDERR from subprocess `{0}`:{1}{2}'.format(
                    p.args, os.linesep, errs.decode('utf-8')
                )
            )
            raise subprocess.CalledProcessError(
                returncode=p.returncode, cmd=p.args, output=outs,
                stderr=errs
            )
