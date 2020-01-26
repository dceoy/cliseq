#!/usr/bin/env python

from pathlib import Path

import luigi

from ..cli.util import print_log
from .base import ShellTask


class TrimAdapters(ShellTask):
    fq_paths = luigi.ListParameter()
    sample_name = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['trim_dir_path']).joinpath(
                        (
                            Path(Path(p).stem).stem
                            if p.endswith(('.fastq.gz', '.fq.gz'))
                            else Path(p).name
                        ) + f'_val_{i + 1}.fq.gz'
                    )
                )
            ) for i, p in enumerate(self.fq_paths)
        ]

    def run(self):
        run_id = self.sample_name
        print_log(f'Trim adapters:\t{run_id}')
        cutadapt = self.cf['cutadapt']
        fastqc = self.cf['fastqc']
        pigz = self.cf['pigz']
        trim_galore = self.cf['trim_galore']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[cutadapt, fastqc, pigz, trim_galore],
            cwd=self.cf['trim_dir_path']
        )
        self.run_shell(
            args=(
                f'set -e && {trim_galore} --cores={n_cpu} --illumina'
                + ''.join([
                    f' {a}' for a in (
                        ['--paired', *self.fq_paths]
                        if len(self.fq_paths) > 1 else self.fq_paths
                    )
                ])
            ),
            input_files=self.fq_paths,
            output_files=[o.path for o in self.output()]
        )


if __name__ == '__main__':
    luigi.run()
