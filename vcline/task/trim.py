#!/usr/bin/env python

from pathlib import Path

import luigi

from ..cli.util import parse_fq_id, print_log
from .base import ShellTask


class TrimAdapters(ShellTask):
    p = luigi.DictParameter()
    priority = 7

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.p['trim_dir_path']).joinpath(
                        (
                            Path(Path(Path(p).name).stem).stem
                            if p.endswith(('.fastq.gz', '.fq.gz'))
                            else Path(p).name
                        ) + f'_val_{i + 1}.fq.gz'
                    )
                )
            ) for i, p in enumerate(self.p['raw_fq_paths'])
        ]

    def run(self):
        fq_id = parse_fq_id(fq_path=self.p['raw_fq_paths'][0])
        print_log(f'Trim adapters:\t{fq_id}')
        cutadapt = self.p['cutadapt']
        fastqc = self.p['fastqc']
        trim_galore = self.p['trim_galore']
        n_cpu = self.p['n_cpu']
        self.bash_c(
            args=[
                f'{cutadapt} --version',
                f'{fastqc} --version',
                f'{trim_galore} --version',
                (
                    f'set -e && {trim_galore} --cores={n_cpu} --illumina '
                    + ('--paired ' if len(self.p['raw_fq_paths']) > 1 else '')
                    + ' '.join(self.p['raw_fq_paths'])
                )
            ],
            input_files=self.p['raw_fq_paths'],
            output_files=[o.path for o in self.output()],
            cwd=self.p['trim_dir_path'], run_id=fq_id,
            log_dir_path=self.p['log_dir_path']
        )


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
