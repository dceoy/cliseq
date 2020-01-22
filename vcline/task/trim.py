#!/usr/bin/env python

from pathlib import Path

import luigi

from ..cli.util import parse_fq_id, print_log
from .base import ShellTask


class TrimAdapters(ShellTask):
    fq_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 7

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['trim_dir_path']).joinpath(
                        (
                            Path(Path(Path(p).name).stem).stem
                            if p.endswith(('.fastq.gz', '.fq.gz'))
                            else Path(p).name
                        ) + f'_val_{i + 1}.fq.gz'
                    )
                )
            ) for i, p in enumerate(self.fq_paths)
        ]

    def run(self):
        run_id = parse_fq_id(fq_path=self.fq_paths[0])
        print_log(f'Trim adapters:\t{run_id}')
        cutadapt = self.cf['cutadapt']
        fastqc = self.cf['fastqc']
        trim_galore = self.cf['trim_galore']
        n_cpu = self.cf['n_cpu_per_worker']
        self.bash_c(
            args=[
                f'{cutadapt} --version',
                f'{fastqc} --version',
                f'{trim_galore} --version',
                (
                    f'set -e && {trim_galore} --cores={n_cpu} --illumina '
                    + ('--paired ' if len(self.fq_paths) > 1 else '')
                    + ' '.join(self.fq_paths)
                )
            ],
            input_files=self.fq_paths,
            output_files=[o.path for o in self.output()],
            cwd=self.cf['trim_dir_path'], run_id=run_id,
            log_dir_path=self.cf['log_dir_path']
        )


if __name__ == '__main__':
    luigi.run()
