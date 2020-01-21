#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .base import ShellTask
from .ref import CreateBWAIndexes, CreateFASTAIndex, parse_ref_id
from .trim import TrimAdapters, parse_fq_id


@requires(CreateBWAIndexes, CreateFASTAIndex, TrimAdapters)
class AlignReads(ShellTask):
    log_dir_path = luigi.Parameter()
    align_dir_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    bwa = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter()
    priority = 10

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__fq_paths = [i.path for i in self.input()]
        self.__fq_id = parse_fq_id(fq_path=self.__fq_paths[0])
        self.__align_cram_path = generate_align_cram_path(
            fq_paths=self.__fq_paths, ref_fa_path=self.ref_fa_path,
            align_dir_path=self.align_dir_path
        )

    def output(self):
        return luigi.LocalTarget(self.__align_cram_path)

    def run(self):
        print_log(f'Align reads:\t{self.__fq_id}')
        self.init_bash(
            run_id=self.__fq_id, run_dir_path=self.align_dir_path,
            log_dir_path=self.log_dir_path
        )
        self.bash_c(
            args=[
                f'{self.bwa} 2>&1 | grep -e "Version:"',
                f'{self.samtools} 2>&1 | grep -e "Version:"',
                (
                    'set -eo pipefail &&'
                    + f' {self.bwa} mem -t {self.n_cpu}'
                    + ' -R @RG\tID:None\tSM:None\tPL:ILLUMINA\tLB:None'
                    + ' '.join([self.ref_fa_path, *self.__fq_paths])
                    + f' | {self.samtools} view -@ {self.n_cpu}'
                    + f' -T {self.ref_fa_path} -CS - -o -'
                    + f' | {self.samtools} sort -@ {self.n_cpu}'
                    + f' -T {self.ref_fa_path} -o {self.__align_cram_path}'
                ),
                (
                    'set -e && '
                    + f'{self.samtools} index -@ {self.n_cpu}'
                    + f' {self.__align_cram_path}'
                )
            ],
            input_files=self.__fq_paths,
            output_files=[(self.__align_cram_path + s) for s in ['', '.crai']]
        )


def generate_align_cram_path(fq_paths, ref_fa_path, align_dir_path):
    return str(
        Path(align_dir_path).joinpath(
            '.'.join([
                parse_fq_id(fq_paths[0]),
                parse_ref_id(ref_fa_path=ref_fa_path), 'cram'
            ])
        )
    )


if __name__ == '__main__':
    luigi.run()
