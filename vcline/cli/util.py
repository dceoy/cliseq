#!/usr/bin/env python

import logging
import os
import re
from pathlib import Path
from urllib.request import urlcleanup, urlopen, urlretrieve

import yaml
from jinja2 import Environment, FileSystemLoader


def print_log(message):
    logger = logging.getLogger(__name__)
    logger.debug(message)
    print(f'>>\t{message}', flush=True)


def is_url(arg):
    return arg.startswith(('https://', 'http://'))


def parse_ref_id(ref_fa_path):
    return Path(Path(ref_fa_path).name).stem


def parse_fq_id(fq_path):
    return re.sub(r'[\._]R[12][\._].*\.[^\.]+$', '', Path(fq_path).name)


def download_file(url, output_path):
    try:
        print_log(f'Retrieve:{os.linesep}{url} => {output_path}')
        urlretrieve(url, filename=output_path)
        urlcleanup()
    except Exception as e:
        if Path(output_path).exists():
            os.remove(output_path)
        raise e


def download_and_merge_files(urls, output_path, mode='wb'):
    try:
        with open(output_path, mode) as f:
            for u in urls:
                print_log(f'Write:{os.linesep}{u}{output_path}')
                f.write(urlopen(u).read())
    except Exception as e:
        if Path(output_path).exists():
            os.remove(output_path)
        raise e


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
