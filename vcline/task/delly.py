#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .base import ShellTask
from .bcftools import bcftools_sort_and_index
from .ref import CreateExclusionIntervalListBED, FetchReferenceFASTA


@requires(PrepareCRAMTumor, PrepareCRAMNormal, FetchReferenceFASTA,
          CreateExclusionIntervalListBED)
class CallStructualVariantsWithDelly(ShellTask):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        run_dir = Path(self.cf['somatic_sv_delly_dir_path']).joinpath(
            create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(run_dir),
            *[
                luigi.LocalTarget(
                    run_dir.joinpath(run_dir.name + f'.delly.{s}')
                ) for s in ['bcf', 'vcf.gz', 'vcf.gz.tbi']
            ]
        ]

    def run(self):
        run_dir = Path(self.output()[0].path)
        run_id = run_dir.name
        self.print_log(f'Call somatic SVs with Delly:\t{run_id}')
        delly = self.cf['delly']
        bcftools = self.cf['bcftools']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_mb = self.cf['memory_mb_per_worker']
        input_cram_paths = [i[0].path for i in self.input()[0:2]]
        fa_path = self.input()[2][0].path
        exclusion_bed_path = self.input()[3][0].path
        output_bcf_path = self.output()[1].path
        output_vcf_path = self.output()[2].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=delly, cwd=run_dir,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet'], env={'OMP_NUM_THREADS': str(n_cpu)}
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
                *input_cram_paths, fa_path, exclusion_bed_path
            ],
            output_files_or_dirs=output_bcf_path
        )
        bcftools_sort_and_index(
            shelltask=self, bcftools=bcftools, input_vcf_path=output_bcf_path,
            output_vcf_path=output_vcf_path, memory_mb=memory_mb, n_cpu=n_cpu
        )


if __name__ == '__main__':
    luigi.run()
