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
            re.sub(r'\.(cr|b)am$', '.\\1am.\\1ai', self.sam_path)
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
                f'set -e && {self.samtools} quickcheck -v {self.sam_path}'
                + f' && {self.samtools} index -@ {self.n_cpu} {self.sam_path}'
            ),
            input_files_or_dirs=self.sam_path,
            output_files_or_dirs=self.output().path
        )


class SamtoolsView(ShellTask):
    input_sam_path = luigi.Parameter()
    output_sam_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    add_args = luigi.Parameter(default='')
    message = luigi.Parameter(default='')
    remove_input = luigi.BoolParameter(default=True)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        return luigi.LocalTarget(self.output_sam_path)

    def run(self):
        run_id = Path(self.input_sam_path).stem
        if self.message:
            self.print_log(f'{self.message}:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=str(Path(self.output_sam_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path), env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -e && '
                + f'{self.samtools} view -@ {self.n_cpu} -T {self.fa_path}'
                + ' -{0}S{1}'.format(
                    ('C' if self.output_sam_path.endswith('.cram') else 'b'),
                    (f' {self.add_args}' if self.add_args else '')
                )
                + f' -o {self.output_sam_path} {self.input_sam_path}'
            ),
            input_files_or_dirs=[
                self.input_sam_path, self.fa_path, f'{self.fa_path}.fai'
            ],
            output_files_or_dirs=self.output_sam_path
        )
        if self.remove_input:
            self.run_shell(
                args=f'rm -f {self.input_sam_path}',
                input_files_or_dirs=[self.input_sam_path, self.output_sam_path]
            )


class ConvertSAMIntoSortedSAM(ShellTask):
    input_sam_path = luigi.Parameter()
    output_sam_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_per_thread = luigi.Parameter(default='768M')
    index_sam = luigi.BoolParameter(default=True)
    remove_input = luigi.BoolParameter(default=True)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        if self.index_sam:
            return [
                luigi.LocalTarget(p) for p in [
                    self.output_sam_path,
                    re.sub(
                        r'\.(cr|b)am$', '.\\1am.\\1ai', self.output_sam_path
                    )
                ]
            ]
        else:
            return [luigi.LocalTarget(self.output_sam_path)]

    def run(self):
        run_id = Path(self.output_sam_path).stem
        self.print_log(f'Sort SAM:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=str(Path(self.output_sam_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path), env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                'set -eo pipefail && '
                + f'{self.samtools} sort -@ {self.n_cpu}'
                + f' -m {self.memory_per_thread} -O bam -l 0'
                + f' -T {self.output_sam_path}.sort {self.input_sam_path}'
                + f' | {self.samtools} view -@ {self.n_cpu} -T {self.fa_path}'
                + ' -{}S'.format(
                    'C' if self.output_sam_path.endswith('.cram') else 'b'
                )
                + f' -o {self.output_sam_path} -'
            ),
            input_files_or_dirs=[
                self.input_sam_path, self.fa_path, f'{self.fa_path}.fai'
            ],
            output_files_or_dirs=self.output_sam_path
        )
        if self.index_sam:
            yield SamtoolsIndex(
                sam_path=self.output_sam_path, samtools=self.samtools,
                n_cpu=self.n_cpu, log_dir_path=self.log_dir_path,
                remove_if_failed=self.remove_if_failed
            )
        if self.remove_input:
            self.run_shell(
                args=f'rm -f {self.input_sam_path}',
                input_files_or_dirs=[self.output_sam_path, self.input_sam_path]
            )


class MergeSAMsIntoSortedSAM(ShellTask):
    input_sam_paths = luigi.ListParameter()
    output_sam_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_per_thread = luigi.Parameter(default='768M')
    index_sam = luigi.BoolParameter(default=True)
    remove_input = luigi.BoolParameter(default=True)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        if self.index_sam:
            return [
                luigi.LocalTarget(p) for p in [
                    self.output_sam_path,
                    re.sub(
                        r'\.(cr|b)am$', '.\\1am.\\1ai', self.output_sam_path
                    )
                ]
            ]
        else:
            return [luigi.LocalTarget(self.output_sam_path)]

    def run(self):
        if len(self.input_sam_paths) == 1:
            yield ConvertSAMIntoSortedSAM(
                input_sam_path=self.input_sam_paths[0],
                output_sam_path=self.output_sam_path, fa_path=self.fa_path,
                samtools=self.samtools, n_cpu=self.n_cpu,
                memory_per_thread=self.memory_per_thread,
                index_sam=self.index_sam, remove_input=self.remove_input,
                log_dir_path=self.log_dir_path,
                remove_if_failed=self.remove_if_failed
            )
        else:
            run_id = Path(self.output_sam_path).stem
            self.print_log(f'Merge SAMs:\t{run_id}')
            self.setup_shell(
                run_id=run_id, log_dir_path=(self.log_dir_path or None),
                commands=self.samtools,
                cwd=str(Path(self.output_sam_path).parent),
                remove_if_failed=self.remove_if_failed,
                quiet=bool(self.log_dir_path), env={'REF_CACHE': '.ref_cache'}
            )
            self.run_shell(
                args=(
                    'set -eo pipefail && '
                    + f'{self.samtools} merge -@ {self.n_cpu} -rh - '
                    + ' '.join(self.input_sam_paths)
                    + f' | {self.samtools} sort -@ {self.n_cpu}'
                    + f' -m {self.memory_per_thread} -O bam -l 0'
                    + f' -T {self.output_sam_path}.sort -'
                    + f' | {self.samtools} view -@ {self.n_cpu}'
                    + f' -T {self.fa_path}'
                    + ' -{}S'.format(
                        'C' if self.output_sam_path.endswith('.cram') else 'b'
                    )
                    + f' -o {self.output_sam_path} -'
                ),
                input_files_or_dirs=[
                    *self.input_sam_paths, self.fa_path, f'{self.fa_path}.fai'
                ],
                output_files_or_dirs=self.output_sam_path
            )
            if self.index_sam:
                yield SamtoolsIndex(
                    sam_path=self.output_sam_path, samtools=self.samtools,
                    n_cpu=self.n_cpu, log_dir_path=self.log_dir_path,
                    remove_if_failed=self.remove_if_failed
                )
            if self.remove_input:
                self.run_shell(
                    args=('rm -f ' + ' '.join(self.input_sam_paths)),
                    input_files_or_dirs=[
                        self.output_sam_path, *self.input_sam_paths
                    ]
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
        run_id = Path(Path(self.tsv_path).stem).stem
        self.print_log(f'Index a TAB-delimited position file:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.tabix, cwd=str(Path(self.tsv_path).parent),
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=f'set -e && {self.tabix} -p {self.preset} {self.tsv_path}',
            input_files_or_dirs=self.tsv_path,
            output_files_or_dirs=self.output().path
        )


if __name__ == '__main__':
    luigi.run()
