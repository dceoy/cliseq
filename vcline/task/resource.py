#!/usr/bin/env python

import re
import sys
from pathlib import Path

import luigi

from .base import ShellTask


class DownloadResourceFiles(ShellTask):
    urls = luigi.ListParameter()
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    wget = luigi.Parameter(default='wget')
    pbzip2 = luigi.Parameter(default='pbzip2')
    bgzip = luigi.Parameter(default='bgzip')

    def output(self):
        for u in self.urls:
            p = str(Path(self.dest_dir_path).resolve().joinpath(Path(u).name))
            if p.endswith(('.gz', '.bz2')):
                yield luigi.LocalTarget(p)
            elif p.endswith('.bgz'):
                yield luigi.LocalTarget(re.sub(r'\.bgz$', '.gz', p))
            elif p.endswith('.vcf'):
                yield luigi.LocalTarget(f'{p}.gz')
            else:
                yield luigi.LocalTarget(f'{p}.bz2')

    def run(self):
        self.print_log(f'Download resource files:\t{self.dest_dir_path}')
        self.setup_shell(
            commands=[self.wget, self.pbzip2, self.bgzip],
            cwd=self.dest_dir_path, quiet=False
        )
        for u, o in zip(self.urls, self.output()):
            p = str(Path(self.dest_dir_path).resolve().joinpath(Path(u).name))
            self.run_shell(
                args=f'set -e && {self.wget} -qSL {u} -O {p}',
                output_files_or_dirs=p
            )
            if p != o.path:
                if p.endswith('.bgz'):
                    cmd = f'mv {p} {o.path}'
                elif p.endswith('.vcf'):
                    cmd = f'{self.bgzip} -@ {self.n_cpu} {p}'
                else:
                    cmd = f'{self.pbzip2} -p{self.n_cpu} {p}'
                self.run_shell(
                    args=f'set -e && {cmd}', input_files_or_dirs=p,
                    output_files_or_dirs=o.path
                )


class DownloadFuncotatorDataSources(ShellTask):
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    gatk = luigi.Parameter(default='gatk')

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
        existing_data_keys = {
            Path(o.stem).stem[-1] for o in Path(self.dest_dir_path).iterdir()
            if (
                o.name.startswith('funcotator_dataSources.')
                and o.name.endswith('.tar.gz')
            )
        }
        if not ({'s', 'g'} <= existing_data_keys):
            self.setup_shell(
                commands=self.gatk, cwd=self.dest_dir_path, quiet=False
            )
            for k in ['germline', 'somatic']:
                if k[0] not in existing_data_keys:
                    self.run_shell(
                        args=(
                            f'set -e && {self.gatk}{gatk_opts}'
                            + ' FuncotatorDataSourceDownloader'
                            + ' --validate-integrity'
                            + f' --{k}'
                        )
                    )


class WriteAfOnlyVCF(ShellTask):
    src_path = luigi.Parameter(default='')
    src_url = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    curl = luigi.Parameter(default='curl')
    bgzip = luigi.Parameter(default='bgzip')

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
        assert bool(self.src_path or self.src_url)
        dest_path = self.output().path
        self.print_log(f'Write AF-only VCF:\t{dest_path}')
        self.setup_shell(
            commands=[
                *(list() if self.src_path else [self.curl]), self.bgzip,
                sys.executable
            ],
            cwd=self.dest_dir_path, quiet=False
        )
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
            luigi.LocalTarget(
                str(Path(self.dest_dir_path).joinpath(n).resolve())
            ) for n in [
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
                f'set -e && '
                + f'{self.gatk}{gatk_opts} BedToIntervalList'
                + f' --INPUT {self.bed_path}'
                + f' --OUTPUT {interval_list_path}'
                + f' --SEQUENCE_DICTIONARY {seq_dict_path}'
            ),
            input_files_or_dirs=[self.bed_path, seq_dict_path],
            output_files_or_dirs=interval_list_path
        )


if __name__ == '__main__':
    luigi.run()
