#!/usr/bin/env python

from itertools import chain
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .base import ShellTask
from .ref import (CreateBWAIndexes, CreateFASTAIndex, FetchGenomeFASTA,
                  parse_ref_id)
from .trim import TrimAdapters, parse_fq_id


@requires(FetchGenomeFASTA, CreateBWAIndexes, CreateFASTAIndex, TrimAdapters)
class AlignReads(ShellTask):
    p = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.p['align_dir_path']).joinpath(
                        '.'.join([
                            parse_fq_id(fq_path=self.p['raw_fq_paths'][0]),
                            'trimmed',
                            parse_ref_id(ref_fa_path=self.input()[0].path),
                            'cram'
                        ]) + s
                    )
                )
            ) for s in ['', '.crai']
        ]

    def run(self):
        fq_id = parse_fq_id(fq_path=self.p['raw_fq_paths'][0]) + '.trimmed'
        print_log(f'Align reads:\t{fq_id}')
        bwa = self.p['bwa']
        samtools = self.p['samtools']
        n_cpu = self.p['n_cpu']
        r = '@RG\tID:None\tSM:None\tPL:ILLUMINA\tLB:None'
        fa = self.input()[0].path
        cram = self.output()[0].path
        self.bash_c(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'{samtools} 2>&1 | grep -e "Version:"',
                (
                    f'set -eo pipefail && {bwa} mem -t {n_cpu} -R {r} {fa} '
                    + ' '.join([i.path for i in self.input()[3]])
                    + f' | {samtools} view -@ {n_cpu} -T {fa} -CS - -o -'
                    + f' | {samtools} sort -@ {n_cpu} -T {fa} -o {cram}'
                ),
                f'set -e && {samtools} index -@ {n_cpu} {cram}'
            ],
            input_files=[i.path for i in chain.from_iterable(self.input())],
            output_files=[o.path for o in self.output()],
            cwd=self.p['align_dir_path'], run_id=fq_id,
            log_dir_path=self.p['log_dir_path']
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
