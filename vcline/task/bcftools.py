#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask


class BcftoolsIndex(ShellTask):
    vcf_path = luigi.Parameter()
    bcftools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    tbi = luigi.BoolParameter(default=True)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    quiet = luigi.BoolParameter(default=False)
    priority = 100

    def output(self):
        return luigi.LocalTarget(
            re.sub(
                r'\.(bcf|vcf.gz)$',
                '.\\1.{}'.format('tbi' if self.tbi else 'csi'), self.vcf_path
            )
        )

    def run(self):
        run_id = re.sub(r'.(bcf|vcf.gz)$', '', Path(self.vcf_path).name)
        self.print_log(f'Index VCF/BCF:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.bcftools, cwd=Path(self.vcf_path).parent,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet
        )
        bcftools_index(
            shelltask=self, bcftools=self.bcftools, vcf_path=self.vcf_path,
            n_cpu=self.n_cpu, tbi=self.tbi
        )


def bcftools_index(shelltask, bcftools, vcf_path, n_cpu, tbi=True):
    shelltask.run_shell(
        args=(
            f'set -e && {bcftools} index --threads {n_cpu}'
            + (' --tbi' if tbi else ' --csi')
            + f' {vcf_path}'
        ),
        input_files_or_dirs=vcf_path,
        output_files_or_dirs=re.sub(
            r'\.(bcf|vcf.gz)$', '.\\1.{}'.format('tbi' if tbi else 'csi'),
            vcf_path
        )
    )


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
        return [
            luigi.LocalTarget(
                Path(self.dest_dir_path).joinpath(
                    re.sub(
                        r'\.vcf$', f'.norm.vcf.gz{s}',
                        Path(self.input_vcf_path).stem
                    )
                )
            ) for s in ['', '.tbi']
        ]

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
                'set -eo pipefail && '
                + f'{self.bcftools} reheader --fai {self.fa_path}.fai'
                + f' --threads {self.n_cpu} {self.input_vcf_path}'
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
        bcftools_index(
            shelltask=self, bcftools=self.bcftools, vcf_path=output_vcf_path,
            n_cpu=self.n_cpu, tbi=True
        )


class ConcatenateVCFsIntoSortedVCF(ShellTask):
    input_vcf_paths = luigi.ListParameter()
    output_vcf_path = luigi.Parameter()
    bcftools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.IntParameter(default=(4 * 1024))
    remove_input = luigi.BoolParameter(default=True)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    quiet = luigi.BoolParameter(default=False)
    priority = 90

    def output(self):
        return [
            luigi.LocalTarget(self.output_vcf_path + s) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-2])
        self.print_log(
            f'Concatenate VCFs:\t{run_id}' if len(self.input_vcf_paths) == 1
            else f'Sort VCF:\t{run_id}'
        )
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.bcftools, cwd=Path(self.output_vcf_path).parent,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet
        )
        self.run_shell(
            args=(
                (
                    'set -e && '
                    + f'{self.bcftools} sort --max-mem {self.memory_mb}M'
                    + f' --temp-dir {output_vcf_path}.sort --output-type z'
                    + f' --output-file {output_vcf_path}'
                    + f' {self.input_vcf_paths[0]}'
                ) if len(self.input_vcf_paths) == 1 else (
                    'set -eo pipefail && '
                    + f'{self.bcftools} concat --threads {self.n_cpu} '
                    + ' '.join(self.input_vcf_paths)
                    + f' | {self.bcftools} sort --max-mem {self.memory_mb}M'
                    + f' --temp-dir {output_vcf_path}.sort --output-type z'
                    + f' --output-file {output_vcf_path} -'
                )
            ),
            input_files_or_dirs=self.input_vcf_paths,
            output_files_or_dirs=output_vcf_path
        )
        bcftools_index(
            shelltask=self, bcftools=self.bcftools, vcf_path=output_vcf_path,
            n_cpu=self.n_cpu, tbi=True
        )
        if self.remove_input:
            self.run_shell(
                args=('rm -f ' + ' '.join(self.input_vcf_paths)),
                input_files_or_dirs=self.input_vcf_paths
            )


def bcftools_concat_and_index(shelltask, bcftools, input_vcf_paths,
                              output_vcf_path, n_cpu=1, memory_mb=1024):
    bcftools_concat(
        shelltask=shelltask, bcftools=bcftools,
        input_vcf_paths=input_vcf_paths, output_vcf_path=output_vcf_path,
        n_cpu=n_cpu, memory_mb=memory_mb
    )
    bcftools_index(
        shelltask=shelltask, bcftools=bcftools, vcf_path=output_vcf_path,
        n_cpu=n_cpu, tbi=True
    )


def bcftools_concat(shelltask, bcftools, input_vcf_paths, output_vcf_path,
                    n_cpu=1, memory_mb=1024):
    shelltask.run_shell(
        args=(
            'set -eo pipefail && '
            + f'{bcftools} concat --threads {n_cpu} '
            + ' '.join(input_vcf_paths)
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
    bcftools_index(
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
