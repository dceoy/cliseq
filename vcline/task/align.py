#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .base import ShellTask
from .ref import CreateBWAIndexes, CreateFASTAIndex, parse_ref_id
from .trim import TrimAdapters, generate_trim_fq_paths, parse_fq_id


@requires(CreateBWAIndexes, CreateFASTAIndex, TrimAdapters)
class AlignReads(ShellTask):
    params = luigi.DictParameter()
    priority = 10

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.trim_fq_paths = generate_trim_fq_paths(
            fq_paths=self.params['raw_fq_paths'],
            trim_dir_path=self.params['trim_dir_path']
        )
        self.ref_fa_path = self.params['ref_fa_path']
        self.align_cram_path = generate_align_cram_path(
            fq_paths=self.trim_fq_paths,
            ref_fa_path=self.params['ref_fa_path'],
            align_dir_path=self.params['align_dir_path']
        )

    def output(self):
        return luigi.LocalTarget([
            (self.align_cram_path + s) for s in ['', '.crai']
        ])

    def run(self):
        fq_id = parse_fq_id(fq_path=self.trim_fq_paths[0])
        bwa = self.params['bwa']
        samtools = self.params['samtools']
        n_cpu = self.params['n_cpu']
        print_log(f'Align reads:\t{fq_id}')
        self.init_bash(
            run_id=fq_id, run_dir_path=self.params['align_dir_path'],
            log_dir_path=self.params['log_dir_path']
        )
        self.bash_c(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'{samtools} 2>&1 | grep -e "Version:"',
                (
                    'set -eo pipefail &&'
                    + f' {bwa} mem -t {n_cpu}'
                    + ' -R @RG\tID:None\tSM:None\tPL:ILLUMINA\tLB:None'
                    + ' '.join([self.ref_fa_path, *self.trim_fq_paths])
                    + f' | {samtools} view -@ {n_cpu}'
                    + f' -T {self.ref_fa_path} -CS - -o -'
                    + f' | {samtools} sort -@ {n_cpu}'
                    + f' -T {self.ref_fa_path} -o {self.align_cram_path}'
                ),
                (
                    'set -e && '
                    + f'{samtools} index -@ {n_cpu} {self.align_cram_path}'
                )
            ],
            input_files=self.trim_fq_paths,
            output_files=[(self.align_cram_path + s) for s in ['', '.crai']]
        )


def generate_align_cram_path(fq_paths, ref_fa_path, align_dir_path):
    return str(
        Path(align_dir_path).joinpath(
            '.'.join([
                parse_fq_id(fq_paths[0]), 'trimmed',
                parse_ref_id(ref_fa_path=ref_fa_path), 'cram'
            ])
        )
    )


if __name__ == '__main__':
    luigi.run()
