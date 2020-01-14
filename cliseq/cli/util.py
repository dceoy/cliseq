#!/usr/bin/env python

import logging
import shutil
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader


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
