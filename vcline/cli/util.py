#!/usr/bin/env python

import logging
import os
import re
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader


def print_log(message):
    logger = logging.getLogger(__name__)
    logger.debug(message)
    print(f'>>\t{message}', flush=True)


def is_url(arg):
    return arg.startswith(('https://', 'http://'))


def parse_fq_id(fq_path):
    return re.sub(r'[\._]R[12][\._].+$', '', Path(fq_path).name)


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
    print_log(f'Render a file:\t{output_path}')
    with open(output_path, 'w') as f:
        f.write(
            Environment(
                loader=FileSystemLoader(
                    str(Path(__file__).parent.joinpath('../template')),
                    encoding='utf8'
                )
            ).get_template(template).render(data) + os.linesep
        )
