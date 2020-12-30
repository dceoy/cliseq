#!/usr/bin/env python

import os
import shutil
from pathlib import Path

from ftarc.cli.util import print_log, read_yml
from jinja2 import Environment, FileSystemLoader


def write_config_yml(path, src_yml='example_vcline.yml'):
    if Path(path).is_file():
        print_log(f'The file exists:\t{path}')
    else:
        print_log(f'Create a config YAML:\t{path}')
        shutil.copyfile(
            str(Path(__file__).parent.joinpath('../static').joinpath(src_yml)),
            Path(path).resolve()
        )


def render_template(template, data, output_path):
    po = (Path(output_path) if isinstance(output_path, str) else output_path)
    print_log(('Overwrite' if po.exists() else 'Render') + f' a file:\t{po}')
    with po.open(mode='w') as f:
        f.write(
            Environment(
                loader=FileSystemLoader(
                    str(Path(__file__).parent.joinpath('../template')),
                    encoding='utf8'
                )
            ).get_template(template).render(data) + os.linesep
        )


def load_default_dict(stem):
    return read_yml(
        path=Path(__file__).parent.parent.joinpath(f'static/{stem}.yml')
    )


def parse_cram_id(cram_path):
    prefix = Path(cram_path).stem
    if '.trim.' not in prefix:
        return prefix
    else:
        t = Path(cram_path).stem.split('.')
        return '.'.join(t[:-(t[::-1].index('trim') + 1)])
