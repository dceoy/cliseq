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
            Path(__file__).parent.parent.joinpath('cli/afonlyvcf.py').resolve()
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
                    f'set -e && {self.gatk}{gatk_opts}'
                    + ' FuncotatorDataSourceDownloader'
                    + f' --{s}'
                    + ' --validate-integrity'
                ) for s in ['germline', 'somatic']
            ]
        )


if __name__ == '__main__':
    luigi.run()
