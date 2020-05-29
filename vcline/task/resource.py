#!/usr/bin/env python

import sys
from pathlib import Path

import luigi

from .base import ShellTask


class WritePassingAfOnlyVCF(ShellTask):
    src_path = luigi.Parameter(default='')
    src_url = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    curl = luigi.Parameter(default='curl')
    bgzip = luigi.Parameter(default='bgzip')

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).joinpath(
                Path(Path(self.src_path or self.src_url).stem).stem
                + '.af-only.vcf.gz'
            )
        )

    def run(self):
        assert bool(self.src_path or self.src_url)
        dest_path = self.output().path
        self.print_log(f'Write passing AF-only VCF:\t{dest_path}')
        self.setup_shell(
            commands=[
                *(list() if self.src_path else [self.curl]), self.bgzip,
                sys.executable
            ],
            cwd=self.dest_dir_path, quiet=False
        )
        pyscript_path = str(
            Path(__file__).resolve().parent.parent.joinpath(
                'script/extract_af_only_vcf.py'
            )
        )
        self.run_shell(
            args=(
                'set -eo pipefail && '
                + (
                    f'{self.bgzip} -@ {self.n_cpu} -dc {self.src_path}'
                    if self.src_path else (
                        f'{self.curl} -LS {self.src_url}'
                        + f' | {self.bgzip} -@ {self.n_cpu} -dc'
                    )
                ) + f' | {sys.executable} {pyscript_path} -'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {dest_path}'
            ),
            input_files_or_dirs=(self.src_path if self.src_path else None),
            output_files_or_dirs=dest_path
        )


class CreateIntervalListWithBED(ShellTask):
    fa_path = luigi.Parameter()
    bed_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    gatk = luigi.Parameter(default='gatk')

    def output(self):
        return [
            luigi.LocalTarget(Path(self.dest_dir_path).joinpath(n)) for n in [
                (Path(self.bed_path).stem + '.interval_list'),
                (Path(self.fa_path).stem + '.dict')
            ]
        ]

    def run(self):
        interval_list_path = self.output()[0].path
        self.print_log(f'Create an interval_list file:\t{interval_list_path}')
        gatk_opts = ' --java-options "{}"'.format(
            ' '.join([
                '-Dsamjdk.compression_level=5',
                '-XX:+UseParallelGC',
                f'-XX:ParallelGCThreads={self.n_cpu}'
            ])
        )
        seq_dict_path = self.output()[1].path
        self.setup_shell(
            commands=self.gatk, cwd=self.dest_dir_path, quiet=False
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{self.gatk}{gatk_opts} CreateSequenceDictionary'
                + f' --REFERENCE {self.fa_path}'
                + f' --OUTPUT {seq_dict_path}'
            ),
            input_files_or_dirs=self.fa_path,
            output_files_or_dirs=seq_dict_path
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk}{gatk_opts} BedToIntervalList'
                + f' --INPUT {self.bed_path}'
                + f' --OUTPUT {interval_list_path}'
                + f' --SEQUENCE_DICTIONARY {seq_dict_path}'
            ),
            input_files_or_dirs=[self.bed_path, seq_dict_path],
            output_files_or_dirs=interval_list_path
        )


if __name__ == '__main__':
    luigi.run()
