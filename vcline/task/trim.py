#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from ..cli.util import print_log
from .base import ShellTask


class TrimAdapters(ShellTask):
    log_dir_path = luigi.Parameter()
    trim_dir_path = luigi.Parameter()
    raw_fq_paths = luigi.ListParameter()
    cutadapt = luigi.Parameter()
    fastqc = luigi.Parameter()
    trim_galore = luigi.Parameter()
    n_cpu = luigi.IntParameter()
    priority = 7

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__fq_id = parse_fq_id(fq_path=self.raw_fq_paths[0])
        self.__trim_fq_paths = generate_trim_fq_paths(
            fq_paths=self.raw_fq_paths, trim_dir_path=self.trim_dir_path
        )

    def output(self):
        return luigi.LocalTarget(self.__trim_fq_paths)

    def run(self):
        print_log(f'Trim adapters:\t{self.__fq_id}')
        self.init_bash(
            run_id=self.__fq_id, run_dir_path=self.trim_dir_path,
            log_dir_path=self.log_dir_path
        )
        fq_args = ' '.join([
            *(['--paired'] if len(self.raw_fq_paths) > 1 else list()),
            *self.raw_fq_paths
        ])
        self.bash_c(
            args=[
                f'{self.cutadapt} --version',
                f'{self.fastqc} --version',
                f'{self.trim_galore} --version',
                f'set -e && {self.trim_galore} --cores={self.n_cpu} {fq_args}'
            ],
            input_files=self.raw_fq_paths, output_files=self.__trim_fq_paths
        )


def parse_fq_id(fq_path):
    return re.sub(r'[\._]R[12][\._].*\.[^\.]+$', '', Path(fq_path).name)


def generate_trim_fq_paths(fq_paths, trim_dir_path):
    return [
        str(
            Path(trim_dir_path).joinpath(
                (
                    Path(Path(Path(p).name).stem).stem
                    if p.endswith(('.fastq.gz', '.fq.gz'))
                    else Path(p).name
                ) + f'_val_{i + 1}.fq.gz'
            )
        ) for i, p in enumerate(fq_paths)
    ]


if __name__ == '__main__':
    luigi.run()
