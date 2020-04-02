#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask


class NormalizeVCF(ShellTask):
    input_vcf_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    bcftools = luigi.Parameter()
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.dest_dir_path).joinpath(
                        re.sub(r'\.vcf$', '', Path(self.input_vcf_path).stem)
                    )
                ) + f'.norm.vcf.gz{s}'
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Normalize VCF:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.bcftools, cwd=self.dest_dir_path,
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=(
                f'set -e && {self.bcftools} norm'
                + f' --fasta-ref {self.fa_path}'
                + ' --check-ref w'
                + ' --rm-dup exact'
                + ' --output-type z'
                + f' --threads {self.n_cpu}'
                + f' --output {output_vcf_path}'
                + f' {self.input_vcf_path}'
            ),
            input_files_or_dirs=[self.input_vcf_path, self.fa_path],
            output_files_or_dirs=output_vcf_path
        )
        self.run_shell(
            args=(
                f'set -e && {self.bcftools} index'
                + f' --tbi --threads {self.n_cpu} {output_vcf_path}'
            ),
            input_files_or_dirs=output_vcf_path,
            output_files_or_dirs=f'{output_vcf_path}.tbi'
        )


if __name__ == '__main__':
    luigi.run()
