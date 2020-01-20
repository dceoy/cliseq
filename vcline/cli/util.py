#!/usr/bin/env python

import logging
import os
from pathlib import Path
from urllib.request import urlcleanup, urlopen, urlretrieve

import yaml


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


def read_yml(path):
    with open(path, 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    return d
