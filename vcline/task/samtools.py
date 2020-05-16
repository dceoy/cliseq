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
    quiet = luigi.BoolParameter(default=False)
    priority = 100

    def output(self):
        return luigi.LocalTarget(self.fa_path + '.fai')

    def run(self):
        run_id = Path(self.fa_path).stem
        self.print_log(f'Index FASTA:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=Path(self.fa_path).parent,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet
        )
        samtools_faidx(
            shelltask=self, samtools=self.samtools, fa_path=self.fa_path
        )


def samtools_faidx(shelltask, samtools, fa_path):
    shelltask.run_shell(
        args=f'set -e && {samtools} faidx {fa_path}',
        input_files_or_dirs=fa_path, output_files_or_dirs=f'{fa_path}.fai'
    )


class SamtoolsIndex(ShellTask):
    sam_path = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    quiet = luigi.BoolParameter(default=False)
    priority = 100

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.(cr|b)am$', '.\\1am.\\1ai', self.sam_path)
        )

    def run(self):
        run_id = Path(self.sam_path).stem
        self.print_log(f'Index CRAM/BAM:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=Path(self.sam_path).parent,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet
        )
        samtools_index(
            shelltask=self, samtools=self.samtools, sam_path=self.sam_path,
            n_cpu=self.n_cpu
        )


def samtools_index(shelltask, samtools, sam_path, n_cpu=1):
    shelltask.run_shell(
        args=(
            f'set -e && {samtools} quickcheck -v {sam_path}'
            + f' && {samtools} index -@ {n_cpu} {sam_path}'
        ),
        input_files_or_dirs=sam_path,
        output_files_or_dirs=re.sub(r'\.(cr|b)am$', '.\\1am.\\1ai', sam_path)
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
    index_sam = luigi.BoolParameter(default=False)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    quiet = luigi.BoolParameter(default=False)
    priority = 100

    def output(self):
        return luigi.LocalTarget(self.output_sam_path)

    def run(self):
        run_id = Path(self.input_sam_path).stem
        if self.message:
            message = self.message
        elif (self.input_sam_path.endswith('.bam')
              and self.output_sam_path.endswith('.cram')):
            message = 'Convert BAM to CRAM'
        elif (self.input_sam_path.endswith('.cram')
              and self.output_sam_path.endswith('.bam')):
            message = 'Convert CRAM to BAM'
        else:
            message = None
        if message:
            self.print_log(f'{message}:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.samtools, cwd=Path(self.output_sam_path).parent,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet,
            env={'REF_CACHE': '.ref_cache'}
        )
        if self.index_sam:
            samtools_view_and_index(
                shelltask=self, samtools=self.samtools,
                input_sam_path=self.input_sam_path, fa_path=self.fa_path,
                output_sam_path=self.output_sam_path, n_cpu=self.n_cpu,
                add_args=self.add_args
            )
        else:
            samtools_view(
                shelltask=self, samtools=self.samtools,
                input_sam_path=self.input_sam_path, fa_path=self.fa_path,
                output_sam_path=self.output_sam_path, n_cpu=self.n_cpu,
                add_args=self.add_args
            )
        if self.remove_input:
            self.run_shell(
                args=f'rm -f {self.input_sam_path}',
                input_files_or_dirs=self.input_sam_path
            )


def samtools_view_and_index(shelltask, samtools, input_sam_path, fa_path,
                            output_sam_path, n_cpu=1, add_args=None):
    samtools_view(
        shelltask=shelltask, samtools=samtools, input_sam_path=input_sam_path,
        fa_path=fa_path, output_sam_path=output_sam_path, n_cpu=n_cpu,
        add_args=add_args
    )
    samtools_index(
        shelltask=shelltask, samtools=samtools, sam_path=output_sam_path,
        n_cpu=n_cpu
    )


def samtools_view(shelltask, samtools, input_sam_path, fa_path,
                  output_sam_path, n_cpu=1, add_args=None):
    shelltask.run_shell(
        args=(
            f'set -e && {samtools} quickcheck -v {input_sam_path}'
            + f' && {samtools} view -@ {n_cpu} -T {fa_path}'
            + ' -{0}S{1}'.format(
                ('C' if output_sam_path.endswith('.cram') else 'b'),
                (f' {add_args}' if add_args else '')
            )
            + f' -o {output_sam_path} {input_sam_path}'
        ),
        input_files_or_dirs=[
            input_sam_path, fa_path, f'{fa_path}.fai'
        ],
        output_files_or_dirs=output_sam_path
    )


class SortSAM(ShellTask):
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
    quiet = luigi.BoolParameter(default=False)
    priority = 90

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
            commands=self.samtools, cwd=Path(self.output_sam_path).parent,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet,
            env={'REF_CACHE': '.ref_cache'}
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
            samtools_index(
                shelltask=self, samtools=self.samtools,
                sam_path=self.output_sam_path, n_cpu=self.n_cpu
            )
        if self.remove_input:
            self.run_shell(
                args=f'rm -f {self.input_sam_path}',
                input_files_or_dirs=self.input_sam_path
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
    quiet = luigi.BoolParameter(default=False)
    priority = 80

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
            yield SortSAM(
                input_sam_path=self.input_sam_paths[0],
                output_sam_path=self.output_sam_path, fa_path=self.fa_path,
                samtools=self.samtools, n_cpu=self.n_cpu,
                memory_per_thread=self.memory_per_thread,
                index_sam=self.index_sam, remove_input=self.remove_input,
                log_dir_path=self.log_dir_path,
                remove_if_failed=self.remove_if_failed, quiet=self.quiet
            )
        else:
            run_id = Path(self.output_sam_path).stem
            self.print_log(f'Merge SAMs:\t{run_id}')
            self.setup_shell(
                run_id=run_id, log_dir_path=(self.log_dir_path or None),
                commands=self.samtools, cwd=Path(self.output_sam_path).parent,
                remove_if_failed=self.remove_if_failed, quiet=self.quiet,
                env={'REF_CACHE': '.ref_cache'}
            )
            self.run_shell(
                args=(
                    'set -eo pipefail && '
                    + f'{self.samtools} merge -@ {self.n_cpu} -r - '
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
                samtools_index(
                    shelltask=self, samtools=self.samtools,
                    sam_path=self.output_sam_path, n_cpu=self.n_cpu
                )
            if self.remove_input:
                self.run_shell(
                    args=('rm -f ' + ' '.join(self.input_sam_paths)),
                    input_files_or_dirs=self.input_sam_paths
                )


def samtools_merge_and_index(shelltask, samtools, input_sam_paths, fa_path,
                             output_sam_path, n_cpu=1, memory_mb=1024):
    samtools_merge(
        shelltask=shelltask, samtools=samtools,
        input_sam_paths=input_sam_paths, fa_path=fa_path,
        output_sam_path=output_sam_path, n_cpu=n_cpu, memory_mb=memory_mb
    )
    samtools_index(
        shelltask=shelltask, samtools=samtools, sam_path=output_sam_path,
        n_cpu=n_cpu
    )


def samtools_merge(shelltask, samtools, input_sam_paths, fa_path,
                   output_sam_path, n_cpu=1, memory_mb=1024):
    memory_mb_per_thread = int(memory_mb / n_cpu / 2)
    shelltask.run_shell(
        args=(
            'set -eo pipefail && {samtools} merge -@ {n_cpu} -r - '
            + ' '.join(input_sam_paths)
            + f' | {samtools} sort -@ {n_cpu} -m {memory_mb_per_thread}M'
            + f' -O bam -l 0 -T {output_sam_path}.sort -'
            + f' | {samtools} view -@ {n_cpu}'
            + f' -T {fa_path}'
            + ' -{}S'.format(
                'C' if output_sam_path.endswith('.cram') else 'b'
            )
            + f' -o {output_sam_path} -'
        ),
        input_files_or_dirs=[
            *input_sam_paths, fa_path, f'{fa_path}.fai'
        ],
        output_files_or_dirs=output_sam_path
    )


if __name__ == '__main__':
    luigi.run()
