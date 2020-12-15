#!/usr/bin/env python

import re
import sys
from pathlib import Path

import luigi
from ftarc.task.base import ShellTask
from ftarc.task.downloader import DownloadResourceFiles


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
        return bool(
            {'s', 'g'} <= set(self._fetch_existing_funcotator_data().keys())
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


class DownloadSnpEffDataSource(ShellTask):
    dest_dir_path = luigi.Parameter(default='.')
    genome_version = luigi.Parameter(default='GRCh38')
    snpeff = luigi.Parameter(default='snpEff')

    def output(self):
        dir_data_paths = self._fetch_existing_snpeff_data()
        if dir_data_paths:
            return luigi.LocalTarget(dir_data_paths[0])
        else:
            return super().output()

    def complete(self):
        return bool(self._fetch_existing_snpeff_data())

    def _fetch_existing_snpeff_data(self):
        dest_dir = Path(self.dest_dir_path)
        snpeff_config = dest_dir.joinpath('snpEff.config')
        if not snpeff_config.is_file():
            return list()
        else:
            data_dir = dest_dir.joinpath('data')
            with open(snpeff_config, 'r') as f:
                for s in f:
                    if re.match(r'\s*data\.dir\s*=\s*\./data/', s):
                        data_dir = dest_dir.joinpath(
                            re.sub(r'\s*data\.dir\s*=\s*', '', s.strip())
                        )
                        break
            if not data_dir.is_dir():
                return list()
            else:
                return [
                    str(o) for o in data_dir.iterdir()
                    if o.name.startswith(self.genome_version) and o.is_dir()
                ]

    def run(self):
        self.print_log(f'Download SnpEff data source:\t{self.dest_dir_path}')
        src_config = Path(self.snpeff).parent.parent.joinpath('snpEff.config')
        dest_config = Path(self.dest_dir_path).joinpath(src_config.name)
        self.setup_shell(
            commands=self.snpeff, cwd=self.dest_dir_path, quiet=False
        )
        self.run_shell(
            args=[
                f'set -e && cp {src_config} {dest_config}',
                (
                    'set -e && '
                    + f'{self.snpeff} databases'
                    + f' | grep -e "^{self.genome_version}[\\.0-9]*\\s"'
                    + ' | cut -f 1'
                    + f' | xargs {self.snpeff} download'
                    + f' -verbose -config {dest_config}'
                )
            ],
            output_files_or_dirs=dest_config
        )


class DownloadAndConvertVCFsIntoPassingAfOnlyVCF(luigi.Task):
    src_url = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    wget = luigi.Parameter(default='wget')
    bgzip = luigi.Parameter(default='bgzip')

    def requires(self):
        return DownloadResourceFiles(
            src_url=[self.src_url], dest_dir_path=self.dest_dir_path,
            n_cpu=self.n_cpu, wget=self.wget, bgzip=self.bgzip
        )

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).joinpath(
                Path(Path(self.input()[0].path).stem).stem + '.af-only.vcf.gz'
            )
        )

    def run(self):
        yield WritePassingAfOnlyVCF(
            src_path=self.input()[0].path, dest_dir_path=self.dest_dir_path,
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
        pyscript = Path(__file__).resolve().parent.parent.joinpath(
            'script/extract_af_only_vcf.py'
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
                ) + f' | {sys.executable} {pyscript} -'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {dest_path}'
            ),
            input_files_or_dirs=(self.src_path if self.src_path else None),
            output_files_or_dirs=dest_path
        )


if __name__ == '__main__':
    luigi.run()
