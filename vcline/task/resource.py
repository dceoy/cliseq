#!/usr/bin/env python

import sys
from pathlib import Path

import luigi

from .base import ShellTask


class WriteAfOnlyVCF(ShellTask):
    src_path = luigi.Parameter(default='')
    src_url = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    curl = luigi.Parameter(default='curl')
    bgzip = luigi.Parameter(default='bgzip')
    n_cpu = luigi.IntParameter(default=1)

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.dest_dir_path).joinpath(
                    Path(Path(self.src_path or self.src_url).stem).stem
                    + '.af-only.vcf.gz'
                ).resolve()
            )
        )

    def run(self):
        assert bool(self.src_url or self.src_path)
        run_id = Path(Path(self.src_url).stem).stem
        self.print_log(f'Write AF-only VCF:\t{run_id}')
        self.setup_shell(
            commands=[self.curl, self.bgzip, sys.executable],
            cwd=self.dest_dir_path, quiet=False
        )
        dest_path = self.output().path
        pyscript_path = str(
            Path(__file__).parent.parent.joinpath(
                'script/extract_af_only_vcf.py'
            ).resolve()
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && '
                + (
                    f'{self.bgzip} -@ {self.n_cpu} -dc {self.src_path}'
                    if self.src_path else (
                        f'{self.curl} -LS {self.src_url}'
                        + f' | {self.bgzip} -@ {self.n_cpu} -dc'
                    )
                ) + f' | {sys.executable} {pyscript_path}'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {dest_path}'
            ),
            input_files=(self.src_path if self.src_path else None),
            output_files=dest_path
        )


class DownloadFuncotatorDataSources(ShellTask):
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)

    def run(self):
        run_id = Path(self.dest_dir_path).name
        self.print_log(f'Download Funcotator data sources:\t{run_id}')
        gatk_opts = ' --java-options "{}"'.format(
            ' '.join([
                '-Dsamjdk.compression_level=5',
                '-XX:+UseParallelGC',
                f'-XX:ParallelGCThreads={self.n_cpu}'
            ])
        )
        self.setup_shell(
            commands=self.gatk, cwd=self.dest_dir_path, quiet=False
        )
        self.run_shell(
            args=[
                (
                    f'set -e && '
                    + f'{self.gatk}{gatk_opts} FuncotatorDataSourceDownloader'
                    + f' --{s}'
                    + ' --validate-integrity'
                ) for s in ['germline', 'somatic']
            ]
        )


class CreateIntervalListWithBED(ShellTask):
    fa_path = luigi.Parameter()
    bed_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)

    def output(self):
        return [
            luigi.LocalTarget(
                str(Path(self.dest_dir_path).joinpath(n).resolve())
            ) for n in [
                (Path(self.bed_path).stem + '.interval_list'),
                (Path(self.fa_path).stem + '.dict')
            ]
        ]

    def run(self):
        interval_list_path = self.output()[0].path
        run_id = Path(interval_list_path).stem
        self.print_log(f'Create an interval_list file:\t{run_id}')
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
            input_files=self.fa_path, output_files=seq_dict_path
        )
        self.run_shell(
            args=(
                f'set -e && '
                + f'{self.gatk}{gatk_opts} BedToIntervalList'
                + f' --INPUT {self.bed_path}'
                + f' --OUTPUT {interval_list_path}'
                + f' --SEQUENCE_DICTIONARY {seq_dict_path}'
            ),
            input_files=[self.bed_path, seq_dict_path],
            output_files=interval_list_path
        )


if __name__ == '__main__':
    luigi.run()
