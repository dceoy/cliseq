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
        self.__sh = None
        self.__cwd = None

    @classmethod
    def setup_bash(cls, run_id=None, log_dir_path=None, work_dir_path=None):
        cls.__sh = ShellOperator(
            log_txt=(
                str(
                    Path(log_dir_path or '.').joinpath(
                        f'{cls.__module__}.{cls.__name__}.{run_id}.sh.log.txt'
                    ).resolve()
                ) if run_id else None
            ),
            quiet=True, clear_log_txt=False,
            logger=logging.getLogger(__name__), print_command=True,
            executable='/bin/bash'
        )
        cls.__cwd = work_dir_path

    @classmethod
    def run_bash(cls, *args, **kwargs):
        if kwargs.get('cwd'):
            cls.__cwd = kwargs['cwd']
        cls.__sh.run(
            *args, **{k: v for k, v in kwargs.items() if k != 'cwd'},
            cwd=cls.__cwd
        )
