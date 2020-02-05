#!/usr/bin/env python

import logging
from datetime import timedelta
from pathlib import Path

import coloredlogs
import luigi
from shoper.shelloperator import ShellOperator


class BaseTask(luigi.Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        coloredlogs.install(level=logging.root.level)

    @luigi.Task.event_handler(luigi.Event.PROCESSING_TIME)
    def print_execution_time(self, processing_time):
        logger = logging.getLogger(__name__)
        message = '{0}.{1}:\t{2} elapsed.'.format(
            self.__class__.__module__, self.__class__.__name__,
            timedelta(seconds=processing_time)
        )
        logger.info(message)
        print(message, flush=True)


class ShellTask(BaseTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__sh = None
        self.__run_kwargs = None

    @classmethod
    def setup_shell(cls, run_id=None, log_dir_path=None, commands=None,
                    clear_log_txt=False, print_command=True, quiet=True,
                    executable='/bin/bash', **run_kwargs):
        cls.__sh = ShellOperator(
            log_txt=(
                str(
                    Path(log_dir_path or '.').joinpath(
                        f'{cls.__module__}.{cls.__name__}.{run_id}.sh.log.txt'
                    ).resolve()
                ) if run_id else None
            ),
            quiet=quiet, clear_log_txt=clear_log_txt,
            logger=logging.getLogger(__name__), print_command=print_command,
            executable=executable
        )
        cls.__run_kwargs = run_kwargs
        if commands:
            cls.run_shell(args=list(cls._generate_version_commands(commands)))

    @classmethod
    def run_shell(cls, *args, **kwargs):
        cls.__sh.run(
            *args, **kwargs,
            **{k: v for k, v in cls.__run_kwargs.items() if k not in kwargs}
        )
        if 'asynchronous' in kwargs:
            cls.__sh.wait()

    @staticmethod
    def _generate_version_commands(commands):
        for c in ([commands] if isinstance(commands, str) else commands):
            if Path(c).name == 'bwa':
                yield f'{c} 2>&1 | grep -e "Program:" -e "Version:"'
            else:
                yield f'{c} --version'
