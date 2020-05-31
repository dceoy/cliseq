#!/usr/bin/env python

import re
import sys
from pathlib import Path

import luigi

from .base import ShellTask


class DownloadResourceFile(ShellTask):
    src_url = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    dest_path = luigi.Parameter(default='')
    n_cpu = luigi.IntParameter(default=1)
    wget = luigi.Parameter(default='wget')
    pbzip2 = luigi.Parameter(default='pbzip2')
    bgzip = luigi.Parameter(default='bgzip')

    def output(self):
        if self.dest_path:
            return luigi.LocalTarget(self.dest_path)
        else:
            p = str(Path(self.dest_dir_path).joinpath(Path(self.src_url).name))
            if p.endswith('.bgz'):
                return luigi.LocalTarget(re.sub(r'\.bgz$', '.gz', p))
            elif p.endswith(('.vcf', '.bed')):
                return luigi.LocalTarget(f'{p}.gz')
            else:
                return luigi.LocalTarget(p)

    def run(self):
        dest_file = Path(self.output().path)
        self.print_log(f'Download resource files:\t{dest_file}')
        if (self.src_url.endswith(('.bgz', '.gz', '.bz2'))
                or (dest_file.suffix == Path(self.src_url).suffix)):
            tmp_path = str(dest_file)
            commands = self.wget
            postproc_args = None
        elif dest_file.name.endswith('.gz'):
            tmp_path = re.sub(r'\.gz$', '', str(dest_file))
            commands = [self.wget, self.bgzip]
            postproc_args = f'{self.bgzip} -@ {self.n_cpu} {tmp_path}'
        elif dest_file.name.endswith('.bz2'):
            tmp_path = re.sub(r'\.bz2$', '', str(dest_file))
            commands = [self.wget, self.pbzip2]
            postproc_args = f'{self.pbzip2} -p{self.n_cpu} {tmp_path}'
        else:
            raise ValueError(f'invalid dest_path: {self.dest_path}')
        self.setup_shell(commands=commands, cwd=dest_file.parent, quiet=False)
        self.run_shell(
            args=f'set -e && {self.wget} -qSL {self.src_url} -O {tmp_path}',
            output_files_or_dirs=tmp_path
        )
        if postproc_args:
            self.run_shell(
                args=f'set -e && {postproc_args}',
                input_files_or_dirs=tmp_path, output_files_or_dirs=dest_file
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
        return [
            str(o) for o in Path(self.dest_dir_path).iterdir()
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
        return DownloadResourceFile(
            src_url=self.src_url, dest_dir_path=self.dest_dir_path,
            n_cpu=self.n_cpu, wget=self.wget, bgzip=self.bgzip
        )

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).joinpath(
                Path(Path(self.input().path).stem).stem + '.af-only.vcf.gz'
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
