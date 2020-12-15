#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from ftarc.task.base import ShellTask
from ftarc.task.samtools import samtools_index


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
    memory_mb_per_thread = int(memory_mb / n_cpu / 8)
    shelltask.run_shell(
        args=(
            f'set -eo pipefail && {samtools} merge -@ {n_cpu} -r -'
            + ''.join([f' {p}' for p in input_sam_paths])
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
