#!/usr/bin/env python

import os
from pathlib import Path

import luigi
from shoper.shelloperator import ShellOperator

from ..cli.util import print_log


class TrimAdapters(luigi.Task):
    run_id = luigi.Parameter()
    foreground_fq = luigi.ListParameter()
    background_fq = luigi.ListParameter()
    trim_dir_path = luigi.Parameter()
    cutadapt = luigi.Parameter()
    fastqc = luigi.Parameter()
    trim_galore = luigi.Parameter()
    n_cpu = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__log_dir = Path(self.log_dir_path)
        self.__trim_dir = Path(self.trim_dir_path)
        self.__trim_fq_paths = {
            k: [
                (
                    (
                        Path(Path(Path(f).name).stem).stem
                        if f.endswith(('.fastq.gz', '.fq.gz')) else f
                    ) + f'_val_{i}.fq.gz'
                ) for i, f in enumerate(getattr(self, k))
            ] for k in ['foreground_fq', 'background_fq']
        }

    def output(self):
        return luigi.LocalTarget(self.__trim_fq_paths)

    def run(self):
        sh = ShellOperator(
            log_txt=str(
                self.__log_dir.joinpath(
                    'trim_adapters.{}.sh.log.txt'.format(self.run_id)
                )
            ),
            quiet=True
        )
        self.__log_dir.mkdir(exist_ok=True)
        self.__trim_dir.mkdir(exist_ok=True)
        for k in ['foreground_fq', 'background_fq']:
            fqs = getattr(self, k)
            fq_args = '--paired {fqs[0]} {fsq[1]}' if len(fqs) > 1 else fqs[0]
            print_log(os.linesep.join(['Trim adapters:'] + fqs))
            sh.run(
                args=[
                    f'{self.cutadapt} --version',
                    f'{self.fastqc} --version',
                    f'{self.trim_galore} --version',
                    f'{self.trim_galore} --cores={self.n_cpu} {fq_args}'
                ],
                input_files=fqs, output_files=self.__trim_fq_paths[k],
                cwd=self.trim_dir_path
            )


if __name__ == '__main__':
    luigi.run()
