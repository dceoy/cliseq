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


def write_config_yml(path):
    if path.is_file():
        print_log('The file exists: {}'.format(path))
    else:
        print_log('Create a config YAML: {}'.format(path))
        shutil.copyfile(
            str(Path(__file__).parent.joinpath('../static/cliseq.yml')),
            Path(path).resolve()
        )


def fetch_executable(cmd):
    executables = [
        cp for cp
        in [Path(p).joinpath(cmd) for p in ':'.split(os.environ['PATH'])]
        if os.access(cp, os.X_OK)
    ]
    if executables:
        return executables[0]
    else:
        raise RuntimeError(f'command not found: {cmd}')


def download_hg(hg_ver='hg38', work_dir='.'):
    print_log('Download genome FASTA: {}'.format(hg_ver))
    urls = read_yml(
        path=str(Path(__file__).parent.joinpath('../static/urls.yml'))
    )['ref_fa_gz'][hg_ver]
    assert all([u.endswith('.gz') for u in urls]), 'invalid gzip URLs'
    output_path = Path(work_dir).joinpath(
        '{}.fa.gz'.format('.'.join([Path(u).stem.stem.name for u in urls]))
    )
    if len(urls) == 1:
        retrieve_url(url=urls[0], output_path=output_path)
    else:
        try:
            with open(output_path, 'wb') as f:
                for u in urls:
                    print_log(
                        'Write:{0}{1} => {2}'.format(
                            os.linesep, u, output_path
                        )
                    )
                    f.write(urlopen(u).read())
        except Exception as e:
            if Path(output_path).exists():
                os.remove(output_path)
            raise e
    print_log('Save: {}'.format(output_path))


def retrieve_url(url, output_path):
    try:
        print_log(
            'Retrieve:{0}{1} => {2}'.format(
                os.linesep, url, output_path
            )
        )
        urlretrieve(url, filename=output_path)
        urlcleanup()
    except Exception as e:
        if Path(output_path).exists():
            os.remove(output_path)
        raise e


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
            print_log('Remove: {}'.format(p))
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


def print_log(message):
    logger = logging.getLogger(__name__)
    logger.info(message)
    print('>>\t{}'.format(message), flush=True)


def set_log_config(debug=None, info=None):
    if debug:
        lv = logging.DEBUG
    elif info:
        lv = logging.INFO
    else:
        lv = logging.WARNING
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S', level=lv
    )
