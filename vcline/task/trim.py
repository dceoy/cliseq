#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask


class PrepareFASTQs(ShellTask):
    fq_paths = luigi.ListParameter()
    sample_name = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        if self.cf['adapter_removal']:
            return [
                luigi.LocalTarget(o) for o in _generate_trimmed_fqs(
                    raw_fq_paths=self.fq_paths,
                    dest_dir_path=self.cf['trim_dir_path']
                )
            ]
        else:
            return [
                luigi.LocalTarget(
                    Path(self.cf['align_dir_path']).joinpath(
                        Path(p).stem + '.gz'
                    ) if p.endswith('.bz2') else p
                ) for p in self.fq_paths
            ]

    def run(self):
        if self.cf['adapter_removal']:
            yield TrimAdapters(
                fq_paths=self.fq_paths, sample_name=self.sample_name,
                cf=self.cf
            )
        else:
            yield [
                Bunzip2AndGzip(bz2_path=p, gz_path=o.path, cf=self.cf)
                for p, o in zip(self.fq_paths, self.output())
                if p.endswith('.bz2')
            ]


class TrimAdapters(ShellTask):
    fq_paths = luigi.ListParameter()
    sample_name = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(o) for o in _generate_trimmed_fqs(
                raw_fq_paths=self.fq_paths,
                dest_dir_path=self.cf['trim_dir_path']
            )
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
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        for i, o in tmp_fq_paths.items():
            self.run_shell(
                args=(
                    'set -eo pipefail && '
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
        if tmp_fq_paths:
            self.run_shell(
                args=('rm -f ' + ' '.join(tmp_fq_paths.values())),
                input_files_or_dirs=list(tmp_fq_paths.values())
            )


class Bunzip2AndGzip(ShellTask):
    bz2_path = luigi.Parameter()
    gz_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(self.gz_path)

    def run(self):
        run_id = Path(self.bz2_path).stem
        self.print_log(f'Bunzip2 and Gzip a file:\t{run_id}')
        pigz = self.cf['pigz']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[pigz, pbzip2], cwd=Path(self.gz_path).parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {pbzip2} -p{n_cpu} -dc {self.bz2_path}'
                + f' | {pigz} -p {n_cpu} -c - > {self.gz_path}'
            ),
            input_files_or_dirs=self.bz2_path,
            output_files_or_dirs=self.gz_path
        )


def _generate_trimmed_fqs(raw_fq_paths, dest_dir_path):
    for i, p in enumerate(raw_fq_paths):
        yield Path(dest_dir_path).joinpath(
            re.sub(
                r'\.(fastq|fq)\.(gz|bz2)$', f'_val_{i + 1}.fq.gz', Path(p).name
            )
        )


if __name__ == '__main__':
    luigi.run()
