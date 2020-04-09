#!/usr/bin/env python

import re
import sys
from pathlib import Path

import luigi

from .base import ShellTask


class DownloadResourceFile(ShellTask):
    src_url = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    wget = luigi.Parameter(default='wget')
    pbzip2 = luigi.Parameter(default='pbzip2')
    bgzip = luigi.Parameter(default='bgzip')

    def output(self):
        p = str(Path(self.dest_dir_path).joinpath(Path(self.src_url).name))
        if p.endswith(('.gz', '.bz2')):
            return luigi.LocalTarget(p)
        elif p.endswith('.bgz'):
            return luigi.LocalTarget(re.sub(r'\.bgz$', '.gz', p))
        elif p.endswith('.vcf'):
            return luigi.LocalTarget(f'{p}.gz')
        else:
            return luigi.LocalTarget(f'{p}.bz2')

    def run(self):
        dest_path = self.output().path
        self.print_log(f'Download resource files:\t{dest_path}')
        if self.src_url.endswith('.bgz'):
            tmp_path = dest_path
            commands = self.wget
            postproc_args = None
        else:
            tmp_path = str(
                Path(dest_path).parent.joinpath(Path(self.src_url).name)
            )
            if dest_path == tmp_path:
                commands = self.wget
                postproc_args = None
            elif tmp_path.endswith('.vcf'):
                commands = [self.wget, self.bgzip]
                postproc_args = f'{self.bgzip} -@ {self.n_cpu} {tmp_path}'
            else:
                commands = [self.wget, self.pbzip2]
                postproc_args = f'{self.pbzip2} -p{self.n_cpu} {tmp_path}'
        self.setup_shell(
            commands=commands, cwd=self.dest_dir_path, quiet=False
        )
        self.run_shell(
            args=f'set -e && {self.wget} -qSL {self.src_url} -O {tmp_path}',
            output_files_or_dirs=tmp_path
        )
        if postproc_args:
            self.run_shell(
                args=f'set -e && {postproc_args}',
                input_files_or_dirs=tmp_path, output_files_or_dirs=dest_path
            )


class DownloadFuncotatorDataSources(ShellTask):
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    gatk = luigi.Parameter(default='gatk')

    def output(self):
        dir_data_dict = self._fetch_existing_funcotator_data()
        if ({'s', 'g'} <= set(dir_data_dict.keys())):
            return [luigi.LocalTarget(dir_data_dict[k]) for k in ['s', 'g']]
        else:
            return super().output()

    def complete(self):
        return (
            ({'s', 'g'} <= set(self._fetch_existing_funcotator_data().keys()))
            or super().complete()
        )

    def _fetch_existing_funcotator_data(self):
        return {
            Path(o.stem).stem[-1]: str(o)
            for o in Path(self.dest_dir_path).iterdir() if (
                o.name.startswith('funcotator_dataSources.')
                and o.name.endswith('.tar.gz')
            )
        }

    def run(self):
        self.print_log(
            f'Download Funcotator data sources:\t{self.dest_dir_path}'
        )
        gatk_opts = ' --java-options "{}"'.format(
            ' '.join([
                '-Dsamjdk.compression_level=5',
                '-XX:+UseParallelGC',
                f'-XX:ParallelGCThreads={self.n_cpu}'
            ])
        )
        dir_data_dict = self._fetch_existing_funcotator_data()
        self.setup_shell(
            commands=self.gatk, cwd=self.dest_dir_path, quiet=False
        )
        for k in ['germline', 'somatic']:
            if k[0] not in dir_data_dict:
                self.run_shell(
                    args=(
                        f'set -e && {self.gatk}{gatk_opts}'
                        + ' FuncotatorDataSourceDownloader'
                        + ' --validate-integrity'
                        + f' --{k}'
                    )
                )


class DownloadAndConvertVCFsIntoPassingAfOnlyVCF(ShellTask):
    src_url = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    wget = luigi.Parameter(default='wget')
    bgzip = luigi.Parameter(default='bgzip')

    def requires(self):
        return DownloadResourceFile(
            src_url=self.src_url, dest_dir_path=self.dest_dir_path,
            n_cpu=self.n_cpu, wget=self.wget, bgzip=self.bgzip
        )

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.dest_dir_path).joinpath(
                    Path(Path(self.input().path).stem).stem
                    + '.af-only.vcf.gz'
                )
            )
        )

    def run(self):
        yield WritePassingAfOnlyVCF(
            src_path=self.input().path, dest_dir_path=self.dest_dir_path,
            n_cpu=self.n_cpu, bgzip=self.bgzip
        )


class WritePassingAfOnlyVCF(ShellTask):
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
                )
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
                str(Path(self.dest_dir_path).joinpath(n))
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
