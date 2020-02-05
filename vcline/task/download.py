#!/usr/bin/env python

from pathlib import Path

import luigi

from ..cli.util import print_log
from .base import ShellTask


class DownloadVCFAndExtractAF(ShellTask):
    vcf_gz_url = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    curl = luigi.Parameter(default='curl')
    sed = luigi.Parameter(default='sed')
    bgzip = luigi.Parameter(default='bgzip')
    n_cpu = luigi.IntParameter(default=1)

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.dest_dir_path).joinpath(
                    Path(Path(self.vcf_gz_url).stem).stem + '.af-only.vcf.gz'
                ).resolve()
            )
        )

    def run(self):
        run_id = Path(Path(self.vcf_gz_url).stem).stem
        print_log(f'Download a gnomAD VCF and extract AF:\t{run_id}')
        af_vcf_gz_path = self.output().path
        sed_arg0 = (
            's/^\\([^#\\t]*\\t[^\\t]*\\t[^\\t]*\\t[^\\t]*\\t[^\\t]*\\t[^\\t]*'
            '\\t[^\\t]*\\t\\)[^\\t]*;*\\(AF=[0-9]*\\.[e0-9+-]*\\);*[^\\t]*/'
            '\\1\\2/'
        )
        sed_arg1 = '/\\tAF=0\\.[e0+-]*$/d'
        self.setup_shell(
            commands=[self.curl, self.sed, self.bgzip], cwd=self.dest_dir_path,
            quiet=False
        )
        self.run_shell(
            args=(
                f'set -e && '
                + f'{self.curl} -LS {self.vcf_gz_url}'
                + f' | {self.bgzip} -@ {self.n_cpu} -dc -'
                + f' | {self.sed} -e \'{sed_arg0}\' -e \'{sed_arg1}\''
                + f' | {self.bgzip} -@ {self.n_cpu} -c < /dev/stdin'
                + f' > {af_vcf_gz_path}'
            ),
            output_files=af_vcf_gz_path
        )


if __name__ == '__main__':
    luigi.run()
