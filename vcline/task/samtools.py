#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask


class SamtoolsFaidx(ShellTask):
    fa_path = luigi.Parameter()
    samtools = luigi.Parameter()
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        return luigi.LocalTarget(self.fa_path + '.fai')

    def run(self):
        run_id = Path(self.fa_path).stem
        self.print_log(f'Index FASTA:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=str(Path(self.fa_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=f'set -e && {self.samtools} faidx {self.fa_path}',
            input_files_or_dirs=self.fa_path,
            output_files_or_dirs=self.output().path
        )


class SamtoolsIndex(ShellTask):
    sam_path = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.(cr|b)am$', '.\1am.\1ai', self.sam_path)
        )

    def run(self):
        run_id = Path(self.sam_path).stem
        self.print_log(f'Index CRAM/BAM:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=str(Path(self.sam_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=(
                f'set -e && {self.samtools} quickcheck -v {self.sam_path} && '
                + f'{self.samtools} index -@ {self.n_cpu} {self.sam_path}'
            ),
            input_files_or_dirs=self.sam_path,
            output_files_or_dirs=self.output().path
        )


class BAM2CRAM(ShellTask):
    bam_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        return luigi.LocalTarget(re.sub(r'\.bam$', '.cram', self.bam_path))

    def run(self):
        run_id = Path(self.bam_path).stem
        self.print_log(f'Convert BAM to CRAM:\t{run_id}')
        cram_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=str(Path(cram_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path), env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -e && {self.samtools} view -@ {self.n_cpu}'
                + f' -T {self.fa_path} -CS -o {cram_path} {self.bam_path}'
            ),
            input_files_or_dirs=[
                self.bam_path, self.fa_path, f'{self.fa_path}.fai'
            ],
            output_files_or_dirs=cram_path
        )
        self.run_shell(
            args=(
                f'set -e && {self.samtools} quickcheck -v {cram_path} && '
                + f'rm -f {self.bam_path}'
            ),
            input_files_or_dirs=[cram_path, self.bam_path]
        )


class MergeBAMsIntoCRAM(ShellTask):
    bam_paths = luigi.ListParameter()
    output_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_per_thread = luigi.Parameter(default='768M')
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        return luigi.LocalTarget(self.output_cram_path)

    def run(self):
        run_id = Path(self.output_cram_path).stem
        self.print_log(f'Merge BAMs:\t{run_id}')
        cram_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=str(Path(cram_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path), env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                'set -eo pipefail && '
                + f'{self.samtools} merge -@ {self.n_cpu} -rh - '
                + ' '.join(self.bam_paths)
                + f' | {self.samtools} sort -@ {self.n_cpu}'
                + f' -m {self.memory_per_thread} -l 0 -T {cram_path}.sort -'
                + f' | {self.samtools} view -@ {self.n_cpu} -T {self.fa_path}'
                + f' -CS -o {cram_path} -'
            ),
            input_files_or_dirs=[
                self.bam_paths, self.fa_path, f'{self.fa_path}.fai'
            ],
            output_files_or_dirs=cram_path
        )
        self.run_shell(
            args=(
                f'set -e && {self.samtools} quickcheck -v {cram_path} && '
                + f'rm -f ' ' '.join(self.bam_paths)
            ),
            input_files_or_dirs=[cram_path, *self.bam_paths]
        )


class Tabix(ShellTask):
    tsv_path = luigi.Parameter()
    tabix = luigi.Parameter()
    preset = luigi.Parameter(default='bed')
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        return luigi.LocalTarget(f'{self.tsv_path}.tbi')

    def run(self):
        run_id = Path(Path(self.sam_path).stem).stem
        self.print_log(f'Index a TAB-delimited position file:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.tabix, cwd=str(Path(self.tsv_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=f'set -e && {self.tabix} -p {self.preset} {self.tsv_path}',
            input_files_or_dirs=self.sam_path,
            output_files_or_dirs=self.output().path
        )


if __name__ == '__main__':
    luigi.run()
