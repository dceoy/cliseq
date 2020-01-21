#!/usr/bin/env python

import logging
from pathlib import Path

import coloredlogs
import luigi
from shoper.shelloperator import ShellOperator


class BaseTask(luigi.Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        coloredlogs.install(
            level=logging.getLevelName(logging.getLogger(__name__))
        )


class ShellTask(BaseTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def init_bash(cls, run_id, run_dir_path, log_dir_path):
        cls.__run_dir = Path(run_dir_path).resolve()
        cls.__log_dir = Path(log_dir_path).resolve()
        cls.sh = ShellOperator(
            log_txt=str(
                cls.__log_dir.joinpath(
                    '.'.join([
                        cls.__module__, cls.__name__, run_id, 'sh.log.txt'
                    ])
                )
            ),
            quiet=True, clear_log_txt=False,
            logger=logging.getLogger(__name__), print_command=True,
            executable='/bin/bash'
        )

    @classmethod
    def bash_c(cls, *args, **kwargs):
        cls.sh.run(*args, **kwargs, cwd=str(cls.__run_dir))
