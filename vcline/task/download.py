#!/usr/bin/env python

from pathlib import Path

import luigi

from ..cli.util import curl_and_write_af_only_vcf_bgz, print_log
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
        self.setup_shell(
            commands=[self.curl, self.bgzip], cwd=self.dest_dir_path,
            quiet=False
        )
        curl_and_write_af_only_vcf_bgz(
            vcf_bgz_url=self.vcf_gz_url, vcf_bgz_path=self.output().path,
            curl=self.curl, bgzip=self.bgzip, n_cpu=self.n_cpu,
            cwd=self.dest_dir_path
        )


if __name__ == '__main__':
    luigi.run()
