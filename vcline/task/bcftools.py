#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from ftarc.task.base import ShellTask


class NormalizeVCF(ShellTask):
    input_vcf_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter()
    bcftools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.IntParameter(default=(4 * 1024))
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    quiet = luigi.BoolParameter(default=False)
    priority = 10

    def output(self):
        output_vcf = Path(self.dest_dir_path).joinpath(
            re.sub(r'\.vcf$', '', Path(self.input_vcf_path).stem)
            + '.norm.vcf.gz'
        )
        return [luigi.LocalTarget(f'{output_vcf}{s}') for s in ['', '.tbi']]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Normalize VCF:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.bcftools, cwd=self.dest_dir_path,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {self.bcftools} reheader'
                + f' --fai {self.fa_path}.fai --threads {self.n_cpu}'
                + f' {self.input_vcf_path}'
                + f' | {self.bcftools} sort --max-mem {self.memory_mb}M'
                + f' --temp-dir {output_vcf_path}.sort -'
                + f' | {self.bcftools} norm --fasta-ref {self.fa_path}'
                + ' --check-ref w --rm-dup exact --output-type z'
                + f' --threads {self.n_cpu} --output {output_vcf_path} -'
            ),
            input_files_or_dirs=[
                self.input_vcf_path, self.fa_path, f'{self.fa_path}.fai'
            ],
            output_files_or_dirs=output_vcf_path
        )
        _bcftools_index(
            shelltask=self, bcftools=self.bcftools, vcf_path=output_vcf_path,
            n_cpu=self.n_cpu, tbi=True
        )


def _bcftools_index(shelltask, bcftools, vcf_path, n_cpu, tbi=True):
    shelltask.run_shell(
        args=(
            f'set -e && {bcftools} index --threads {n_cpu}'
            + (' --tbi' if tbi else ' --csi')
            + f' {vcf_path}'
        ),
        input_files_or_dirs=vcf_path,
        output_files_or_dirs=re.sub(
            r'\.(bcf|vcf.gz)$', '.\\1.{}'.format('tbi' if tbi else 'csi'),
            str(vcf_path)
        )
    )


def bcftools_concat_and_index(shelltask, bcftools, input_vcf_paths,
                              output_vcf_path, n_cpu=1, memory_mb=1024):
    bcftools_concat(
        shelltask=shelltask, bcftools=bcftools,
        input_vcf_paths=input_vcf_paths, output_vcf_path=output_vcf_path,
        n_cpu=n_cpu, memory_mb=memory_mb
    )
    _bcftools_index(
        shelltask=shelltask, bcftools=bcftools, vcf_path=output_vcf_path,
        n_cpu=n_cpu, tbi=True
    )


def bcftools_concat(shelltask, bcftools, input_vcf_paths, output_vcf_path,
                    n_cpu=1, memory_mb=1024):
    shelltask.run_shell(
        args=(
            f'set -eo pipefail && {bcftools} concat --threads {n_cpu}'
            + ''.join([f' {p}' for p in input_vcf_paths])
            + f' | {bcftools} sort --max-mem {memory_mb}M'
            + f' --temp-dir {output_vcf_path}.sort --output-type z'
            + f' --output-file {output_vcf_path} -'
        ),
        input_files_or_dirs=input_vcf_paths,
        output_files_or_dirs=output_vcf_path
    )


def bcftools_sort_and_index(shelltask, bcftools, input_vcf_path,
                            output_vcf_path, memory_mb=1024, n_cpu=1):
    bcftools_sort(
        shelltask=shelltask, bcftools=bcftools, input_vcf_path=input_vcf_path,
        output_vcf_path=output_vcf_path, memory_mb=memory_mb
    )
    _bcftools_index(
        shelltask=shelltask, bcftools=bcftools, vcf_path=output_vcf_path,
        n_cpu=n_cpu, tbi=True
    )


def bcftools_sort(shelltask, bcftools, input_vcf_path, output_vcf_path,
                  memory_mb=1024):
    shelltask.run_shell(
        args=(
            f'set -e && {bcftools} sort --max-mem {memory_mb}M'
            + f' --temp-dir {output_vcf_path}.sort'
            + f' --output-type z --output-file {output_vcf_path}'
            + f' {input_vcf_path}'
        ),
        input_files_or_dirs=input_vcf_path,
        output_files_or_dirs=output_vcf_path
    )


if __name__ == '__main__':
    luigi.run()
