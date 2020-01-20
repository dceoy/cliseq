#!/usr/bin/env python

from itertools import chain
from pathlib import Path

import luigi

from ..cli.util import print_log
from .bash import BashTask


class TrimAdapters(BashTask):
    log_dir_path = luigi.Parameter()
    trim_dir_path = luigi.Parameter()
    run_id = luigi.Parameter()
    fq_paths = luigi.DictParameter()
    cutadapt = luigi.Parameter()
    fastqc = luigi.Parameter()
    trim_galore = luigi.Parameter()
    n_cpu = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__trim_fq_paths = generate_trim_fq_paths(
            fq_paths=self.fq_paths, trim_dir_path=self.trim_dir_path
        )

    def output(self):
        return luigi.LocalTarget(
            list(chain.from_iterable(self.__trim_fq_paths.values()))
        )

    def run(self):
        print_log(f'Trim adapters:\t{self.run_id}')
        self.init_bash(
            log_name=self.run_id, run_dir_path=self.trim_dir_path,
            log_dir_path=self.log_dir_path
        )
        self.bash_c(
            args=[
                f'{self.cutadapt} --version',
                f'{self.fastqc} --version',
                f'{self.trim_galore} --version'
            ]
        )
        for k, v in self.fq_paths.items():
            self.bash_c(
                args=(
                    f'set -e && {self.trim_galore} --cores={self.n_cpu} '
                    + ' '.join(['--paired', *v] if len(v) > 1 else v)
                ),
                input_files=v, output_files=self.__trim_fq_paths[k]
            )


def generate_trim_fq_paths(fq_paths, trim_dir_path):
    return {
        k: [
            str(
                Path(trim_dir_path).joinpath(
                    (
                        Path(Path(Path(p).name).stem).stem
                        if p.endswith(('.fastq.gz', '.fq.gz'))
                        else Path(p).name
                    ) + f'_val_{i + 1}.fq.gz'
                )
            ) for i, p in enumerate(v)
        ] for k, v in fq_paths.items()
    }


if __name__ == '__main__':
    luigi.run()
