#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .bash import BashTask
from .ref import CreateBWAIndexes, CreateFASTAIndex, generate_ref_fa_path
from .trim import TrimAdapters, generate_trim_fq_paths


@requires(CreateBWAIndexes, CreateFASTAIndex, TrimAdapters)
class AlignReads(BashTask):
    log_dir_path = luigi.Parameter()
    ref_dir_path = luigi.Parameter()
    trim_dir_path = luigi.Parameter()
    align_dir_path = luigi.Parameter()
    run_id = luigi.Parameter()
    ref_id = luigi.Parameter()
    fq_paths = luigi.DictParameter()
    bwa = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__trim_fq_paths = generate_trim_fq_paths(
            fq_paths=self.fq_paths, trim_dir_path=self.trim_dir_path
        )
        self.__ref_fa_path = generate_ref_fa_path(
            ref_id=self.ref_id, ref_dir_path=self.ref_dir_path
        )
        self.__align_cram_paths = generate_align_cram_paths(
            trim_fq_paths=self.__trim_fq_paths, ref_id=self.ref_id,
            align_dir_path=self.align_dir_path
        )

    def output(self):
        return luigi.LocalTarget(list(self.__align_cram_paths.values()))

    def run(self):
        print_log(f'Align reads:\t{self.run_id}')
        self.init_bash(
            log_name=self.run_id, run_dir_path=self.align_dir_path,
            log_dir_path=self.log_dir_path
        )
        self.bash_c(
            args=[
                f'{self.bwa} 2>&1 | grep -e "Version:"',
                f'{self.samtools} 2>&1 | grep -e "Version:"'
            ]
        )
        pl = 'ILLUMINA'
        for k, v in self.self.__trim_fq_paths.items():
            cram = self.__align_cram_paths[k]
            self.bash_c(
                args=[
                    (
                        'set -eo pipefail &&'
                        + f' {self.bwa} mem -t {self.n_cpu}'
                        + f' -R @RG\tID:None\tSM:None\tPL:{pl}\tLB:None'
                        + f' {self.__ref_fa_path}'
                        + ' '.join(self.__trim_fq_paths[k])
                        + f' | {self.samtools} view -@ {self.n_cpu}'
                        + f' -T {self.__ref_fa_path} -CS - -o -'
                        + f' | {self.samtools} sort -@ {self.n_cpu}'
                        + f' -T {self.__ref_fa_path} -o {cram}'
                    ),
                    f'set -e && {self.samtools} index -@ {self.n_cpu} {cram}'
                ], input_files=v, output_files=[cram, f'{cram}.crai']
            )


def generate_align_cram_paths(trim_fq_paths, ref_id, align_dir_path):
    return {
        k: str(
            Path(align_dir_path).joinpath(
                re.sub(
                    r'[\._]R[12][\._].*$', '',
                    Path(Path(Path(v[0]).name).stem).stem
                ) + f'.{ref_id}.cram'
            )
        ) for k, v in trim_fq_paths.items()
    }


if __name__ == '__main__':
    luigi.run()
