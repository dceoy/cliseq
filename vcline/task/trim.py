#!/usr/bin/env python

import re
from pathlib import Path

import luigi

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
                        re.sub(r'\.(fastq|fq)\.(gz|bz2)$', '', Path(p).name)
                        + f'_val_{i + 1}.fq.gz'
                    )
                )
            ) for i, p in enumerate(self.fq_paths)
        ]

    def run(self):
        run_id = self.sample_name
        self.print_log(f'Trim adapters:\t{run_id}')
        cutadapt = self.cf['cutadapt']
        fastqc = self.cf['fastqc']
        pigz = self.cf['pigz']
        trim_galore = self.cf['trim_galore']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        tmp_fq_paths = {
            p: '{}.gz'.format(
                Path(self.cf['trim_dir_path']).joinpath(Path(p).stem)
            ) for p in self.fq_paths if p.endswith('.bz2')
        }
        work_fq_paths = [(tmp_fq_paths.get(p) or p) for p in self.fq_paths]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[cutadapt, fastqc, pigz, trim_galore, pbzip2],
            cwd=self.cf['trim_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        for i, o in tmp_fq_paths.items():
            self.run_shell(
                args=(
                    f'set -eo pipefail && '
                    + f'{pbzip2} -p{n_cpu} -dc {i}'
                    + f' | {pigz} -p {n_cpu} -c - > {o}'
                ),
                input_files_or_dirs=i, output_files_or_dirs=o
            )
        self.run_shell(
            args=(
                f'set -e && {trim_galore} --cores={n_cpu} --illumina'
                + ' --output_dir {}'.format(self.cf['trim_dir_path'])
                + ' --fastqc'
                + (' --paired ' if len(work_fq_paths) > 1 else ' ')
                + ' '.join(work_fq_paths)
            ),
            input_files_or_dirs=work_fq_paths,
            output_files_or_dirs=[o.path for o in self.output()]
        )
        if tmp_fq_paths and self.cf['remove_if_failed']:
            self.run_shell(
                args=('rm -f ' + ' '.join(tmp_fq_paths.values())),
                input_files_or_dirs=list(tmp_fq_paths.values())
            )


if __name__ == '__main__':
    luigi.run()
