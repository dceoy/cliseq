#!/usr/bin/env python

import logging
import os
from pathlib import Path
from pprint import pformat

from ..task.ref import CreateBWAIndexes, CreateFASTAIndex
from .util import fetch_executable, print_log, read_yml


def run_analytical_pipeline(config_yml_path, ref_dir_path=None,
                            work_dir_path='.'):
    logger = logging.getLogger(__name__)
    print_log('Run the analytical pipeline: {}'.format(work_dir_path))
    work_dir = Path(work_dir_path)
    assert work_dir.is_dir()

    cf = read_yml(path=Path(config_yml_path).resolve())
    logger.debug('cf:' + os.linesep + pformat(cf))
    assert 'ref_fa' in cf
    assert set(cf['ref_fa']).intersection({'urls', 'paths'})

    cmd_paths = {
        c: fetch_executable(c) for c in
        ['samtools', 'bwa', 'fastqc', 'cutadapt', 'trim_galore']
    }
    ref_dir_path = (
        Path(ref_dir_path).resolve()
        if ref_dir_path else work_dir.joinpath('ref').resolve()
    )
    param_dict = {
        **config_dict, **cmd_paths, 'ref_dir_path': ref_dir_path
    }
    luigi.build(
        [
            CreateFASTAIndex(param_dict=param_dict),
            CreateBWAIndexes(param_dict=param_dict)
        ],
        workers=2, local_scheduler=True
    )
