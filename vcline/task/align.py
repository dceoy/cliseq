#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .base import ShellTask
from .ref import CreateBWAIndexes, CreateFASTAIndex, FetchGenomeFASTA
from .trim import TrimAdapters, parse_fq_id


@requires(FetchGenomeFASTA, CreateBWAIndexes, CreateFASTAIndex, TrimAdapters)
class AlignReads(ShellTask):
    p = luigi.DictParameter()
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.p['align_dir_path']).joinpath(
                    '{0}.trimmed.{1}.cram'.format(
                        parse_fq_id(fq_path=self.input()[3][0].path),
                        Path(Path(self.input()[0].path).name).stem
                    )
                )
            )
        )

    def run(self):
        run_id = Path(Path(self.output().path).name).stem
        print_log(f'Align reads:\t{run_id}')
        bwa = self.p['bwa']
        samtools = self.p['samtools']
        n_cpu = self.p['n_cpu_per_worker']
        r = '\'@RG\\tID:None\\tSM:None\\tPL:ILLUMINA\\tLB:None\''
        fa_path = self.input()[0].path
        index_paths = [o.path for o in [*self.input()[1], self.input()[2]]]
        fq_paths = [i.path for i in self.input()[3]]
        cram_path = self.output().path
        self.bash_c(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'{samtools} 2>&1 | grep -e "Version:"',
                (
                    'set -eo pipefail && '
                    + f'{bwa} mem -t {n_cpu} -R {r} {fa_path} '
                    + ' '.join(fq_paths)
                    + f' | {samtools} view -@ {n_cpu} -T {fa_path} -CS -'
                    + f' -o {cram_path}'
                ),
            ],
            input_files=[fa_path, *index_paths, *fq_paths],
            output_files=cram_path, cwd=self.p['align_dir_path'],
            run_id=run_id, log_dir_path=self.p['log_dir_path']
        )


if __name__ == '__main__':
    luigi.run()
