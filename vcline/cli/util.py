#!/usr/bin/env python

import bz2
import gzip
import logging
import os
import shutil
from pathlib import Path
from urllib.request import urlcleanup, urlopen, urlretrieve

import yaml
from jinja2 import Environment, FileSystemLoader


def print_log(message):
    logger = logging.getLogger(__name__)
    logger.info(message)
    print('>>\t{}'.format(message), flush=True)


def download_file(url, output_path):
    try:
        print_log(
            'Retrieve:{0}{1} => {2}'.format(os.linesep, url, output_path)
        )
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
                print_log(
                    'Write:{0}{1} => {2}'.format(os.linesep, u, output_path)
                )
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


def open_readable_file(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    elif path.endswith('.bz2'):
        return bz2.open(path, 'rt')
    else:
        return open(path, 'r')


def remove_files_if_they_exists(*args):
    for p in args:
        if Path(p).exists():
            print_log('Remove:\t{}'.format(p))
            os.remove(p)


def render_template(output_path, data=None, template=None):
    logger = logging.getLogger(__name__)
    if Path(output_path).is_file():
        logger.info('Skip rendering a file:\t{}'.format(output_path))
    else:
        print_log('Render a file:   \t{}'.format(output_path))
        if data is not None:
            with open(output_path, 'w') as f:
                f.write(
                    Environment(
                        loader=FileSystemLoader(
                            str(Path(__file__).parent.joinpath('../template')),
                            encoding='utf8'
                        )
                    ).get_template(
                        template or (Path(output_path).name + '.j2')
                    ).render(data)
                )
        else:
            shutil.copyfile(
                str(
                    Path(__file__).parent.joinpath('../static').joinpath(
                        template or Path(output_path).name
                    )
                ),
                output_path
            )


def read_yml(path):
    with open(path, 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    return d
