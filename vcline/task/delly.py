#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMs
from .base import ShellTask
from .ref import (CreateExclusionIntervalListBED, CreateFASTAIndex,
                  FetchReferenceFASTA)


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          CreateExclusionIntervalListBED)
class CallStructualVariantsWithDelly(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['delly_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0]]
                        ) + f'.delly.{s}'
                    )
                )
            ) for s in ['bcf', 'vcf.gz', 'vcf.gz.tbi']
        ]

    def run(self):
        output_bcf_path = self.output()[0].path
        run_id = Path(Path(output_bcf_path).stem).stem
        self.print_log(f'Call somatic SVs with Delly:\t{run_id}')
        delly = self.cf['delly']
        bcftools = self.cf['bcftools']
        n_cpu = self.cf['n_cpu_per_worker']
        input_cram_paths = [i[0].path for i in self.input()[0]]
        fa_path = self.input()[1].path
        fai_path = self.input()[2].path
        exclusion_bed_path = self.input()[3][0].path
        output_vcf_path = self.output()[1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[delly, bcftools], cwd=self.cf['delly_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            env={'OMP_NUM_THREADS': str(n_cpu)}
        )
        self.run_shell(
            args=(
                f'set -e && {delly} call'
                + f' --outfile {output_bcf_path}'
                + f' --genome {fa_path}'
                + f' --exclude {exclusion_bed_path}'
                + ''.join([f' {p}' for p in input_cram_paths])
            ),
            input_files_or_dirs=[
                *input_cram_paths, fa_path, fai_path, exclusion_bed_path
            ],
            output_files_or_dirs=output_bcf_path
        )
        self.run_shell(
            args=(
                f'set -e && {bcftools} view'
                + ' --output-type z'
                + f' --output-file {output_vcf_path}'
                + f' --threads {n_cpu}'
                + f' {output_bcf_path}'
            ),
            input_files_or_dirs=output_bcf_path,
            output_files_or_dirs=output_vcf_path
        )
        self.run_shell(
            args=(
                f'set -e && {bcftools} index'
                + ' --tbi'
                + f' --threads {n_cpu}'
                + f' {output_vcf_path}'
            ),
            input_files_or_dirs=output_vcf_path,
            output_files_or_dirs=f'{output_vcf_path}.tbi'
        )


if __name__ == '__main__':
    luigi.run()
