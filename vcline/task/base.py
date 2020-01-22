#!/usr/bin/env python

import logging
from pathlib import Path

import coloredlogs
import luigi
from shoper.shelloperator import ShellOperator


class BaseTask(luigi.Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        coloredlogs.install(level=logging.getLevelName(logging.root.level))


class ShellTask(BaseTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def bash_c(cls, run_id, log_dir_path, *args, **kwargs):
        sh = ShellOperator(
            log_txt=str(
                Path(log_dir_path).joinpath(
                    f'{cls.__module__}.{cls.__name__}.{run_id}.sh.log.txt'
                )
            ),
            quiet=True, clear_log_txt=False,
            logger=logging.getLogger(__name__), print_command=True,
            executable='/bin/bash'
        )
        sh.run(*args, **kwargs)
