#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask


class NormalizeVCF(ShellTask):
    src_vcf_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                re.sub(r'\.vcf\.gz$', '.norm.vcf.gz', self.src_vcf_path)
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Normalize VCF:\t{run_id}')
        bcftools = self.cf['bcftools']
        tabix = self.cf['tabix']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bcftools, cwd=self.cf['bcftools_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {bcftools} norm'
                + f' --fasta-ref {self.fa_path}'
                + ' --check-ref w'
                + ' --rm-dup exact'
                + ' --multiallelics -'
                + ' --output-type z'
                + f' --threads {n_cpu}'
                + f' --output {output_vcf_path}'
                + f' {self.src_vcf_path}'
            ),
            input_files_or_dirs=[self.src_vcf_path, self.fa_path],
            output_files_or_dirs=output_vcf_path
        )
        self.run_shell(
            args=f'set -e && {tabix} -p vcf {output_vcf_path}',
            input_files_or_dirs=output_vcf_path,
            output_files_or_dirs=f'{output_vcf_path}.tbi'
        )


if __name__ == '__main__':
    luigi.run()
