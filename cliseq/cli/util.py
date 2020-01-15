#!/usr/bin/env python

import logging
import os
import shutil
from pathlib import Path
from urllib.error import URLError
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


def download_hg(ver='hg38', work_dir='.'):
    print_log('Download genome FASTA: {}'.format(ver))
    urls = read_yml(
        path=str(Path(__file__).parent.joinpath('../static/urls.yml'))
    )['ref_fa_gz'][ver]
    assert all([u.endswith('.gz') for u in urls]), 'invalid gzip URLs'
    output_path = Path(work_dir).joinpath(
        '{}.fa.gz'.format('.'.join([Path(u).suffix.suffix.name for u in urls]))
    )
    try:
        if len(urls) == 1:
            print_log(
                'Retrieve:{0}{1} => {2}'.format(
                    os.linesep, urls[0], output_path
                )
            )
            urlretrieve(urls[0], filename=output_path)
            urlcleanup()
        else:
            with open(output_path, 'wb') as f:
                for u in urls:
                    print_log(
                        'Write:{0}{1} => {2}'.format(
                            os.linesep, u, output_path
                        )
                    )
                    f.write(urlopen(u).read())
    except URLError as e:
        logger = logging.getLogger(__name__)
        logger.error(e)
        os.remove(output_path)
        raise e
    else:
        print_log('Save: {}'.format(output_path))


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
