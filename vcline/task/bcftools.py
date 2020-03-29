#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from .base import ShellTask
from .ref import FetchReferenceFASTA


@requires(FetchReferenceFASTA)
class NormalizeVCF(ShellTask):
    input_vcf_path = luigi.Parameter()
    output_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(self.output_vcf_path + s) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Normalize VCF:\t{run_id}')
        bcftools = self.cf['bcftools']
        n_cpu = self.cf['n_cpu_per_worker']
        fa_path = self.input().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bcftools, cwd=str(Path(self.output_vcf_path).parent),
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {bcftools} norm'
                + f' --fasta-ref {fa_path}'
                + ' --check-ref w'
                + ' --rm-dup exact'
                + ' --multiallelics -'
                + ' --output-type z'
                + f' --threads {n_cpu}'
                + f' --output-file {output_vcf_path}'
                + f' {self.input_vcf_path}'
            ),
            input_files_or_dirs=[self.input_vcf_path, fa_path],
            output_files_or_dirs=output_vcf_path
        )
        self.run_shell(
            args=(
                f'set -e && {bcftools} index'
                + ' --tbi'
                + f' --threads {n_cpu}'
                + f' {output_vcf_path}'
            ),
            input_files_or_dirs=output_vcf_path,
            output_files_or_dirs=f'{output_vcf_path}.tbi'
        )


if __name__ == '__main__':
    luigi.run()
