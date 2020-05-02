#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask


class NormalizeVCF(ShellTask):
    input_vcf_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter()
    bcftools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.IntParameter(default=(4 * 1024))
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 50

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
                f'set -eo pipefail && '
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
        yield BcftoolsIndex(
            vcf_path=output_vcf_path, bcftools=self.bcftools, n_cpu=self.n_cpu,
            tbi=True, log_dir_path=self.log_dir_path,
            remove_if_failed=self.remove_if_failed
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
    priority = 50

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
            commands=self.bcftools, cwd=str(Path(self.output_vcf_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=(
                (
                    f'set -e && '
                    + f'{self.bcftools} sort --max-mem {self.memory_mb}M'
                    + f' --temp-dir {output_vcf_path}.sort --output-type z'
                    + f' --output-file {output_vcf_path}'
                    + f' {self.input_vcf_paths[0]}'
                ) if len(self.input_vcf_paths) == 1 else (
                    f'set -eo pipefail && '
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
        yield BcftoolsIndex(
            vcf_path=output_vcf_path, bcftools=self.bcftools, n_cpu=self.n_cpu,
            tbi=True, log_dir_path=self.log_dir_path,
            remove_if_failed=self.remove_if_failed
        )
        if self.remove_input:
            self.run_shell(
                args=('rm -f ' + ' '.join(self.input_vcf_paths)),
                input_files_or_dirs=[output_vcf_path, *self.input_vcf_paths],
            )


class BcftoolsIndex(ShellTask):
    vcf_path = luigi.Parameter()
    bcftools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    tbi = luigi.BoolParameter(default=True)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 50

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
            commands=self.bcftools, cwd=str(Path(self.vcf_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=(
                f'set -e && {self.bcftools} index --threads {self.n_cpu}'
                + (' --tbi' if self.tbi else ' --csi')
                + f' {self.vcf_path}'
            ),
            input_files_or_dirs=self.vcf_path,
            output_files_or_dirs=self.output().path
        )


class SortVCF(ShellTask):
    input_vcf_path = luigi.Parameter()
    output_vcf_path = luigi.Parameter()
    bcftools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.IntParameter(default=(4 * 1024))
    index_vcf = luigi.BoolParameter(default=True)
    remove_input = luigi.BoolParameter(default=True)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(f'{self.output_vcf_path}{s}')
            for s in ['', '.tbi']
        ]

    def run(self):
        run_id = Path(Path(self.output_vcf_path).stem).stem
        self.print_log(f'Sort VCF:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.bcftools,
            cwd=str(Path(self.output_vcf_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && '
                + f'{self.bcftools} sort --max-mem {self.memory_mb}M'
                + f' --temp-dir {self.output_vcf_path}.sort'
                + f' --output-type z --output-file {self.output_vcf_path}'
                + f' {self.input_vcf_path}'
            ),
            input_files_or_dirs=self.input_vcf_path,
            output_files_or_dirs=self.output_vcf_path
        )
        if self.index_vcf:
            yield BcftoolsIndex(
                vcf_path=self.output_vcf_path, bcftools=self.bcftools,
                n_cpu=self.n_cpu, tbi=True, log_dir_path=self.log_dir_path,
                remove_if_failed=self.remove_if_failed
            )
        if self.remove_input:
            self.run_shell(
                args=f'rm -f {self.input_vcf_path}',
                input_files_or_dirs=[self.output_vcf_path, self.input_vcf_path]
            )


if __name__ == '__main__':
    luigi.run()
