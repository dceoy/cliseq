#!/usr/bin/env python

import logging
import os
from pathlib import Path
from pprint import pformat

from .util import print_log, read_yml


def run_analytical_pipeline(config_yml_path, ref_dir_path=None, work_dir='.'):
    logger = logging.getLogger(__name__)
    print_log('Run the analytical pipeline: {}'.format(work_dir))
    cf = read_yml(path=Path(config_yml_path).resolve())
    logger.debug('args:' + os.linesep + pformat(cf))
