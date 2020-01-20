#!/usr/bin/env python

import logging
from pathlib import Path

import luigi
from shoper.shelloperator import ShellOperator


class BashTask(luigi.Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def init_bash(cls, log_name, run_dir_path, log_dir_path):
        cls.__run_dir = Path(run_dir_path).resolve()
        log_dir = Path(log_dir_path).resolve()
        cls.sh = ShellOperator(
            log_txt=str(
                log_dir.joinpath(
                    '.'.join([
                        cls.__module__, cls.__name__, log_name, 'sh.log.txt'
                    ])
                )
            ),
            quiet=True, clear_log_txt=False,
            logger=logging.getLogger(__name__), print_command=True,
            executable='/bin/bash'
        )
        log_dir.mkdir(exist_ok=True)

    @classmethod
    def bash_c(cls, *args, **kwargs):
        cls.__run_dir.mkdir(exist_ok=True)
        cls.sh.run(
            *args, **kwargs, cwd=str(cls.__run_dir), prompt=None,
            in_background=False, remove_if_failed=True, remove_previous=False,
            skip_if_exist=True
        )
