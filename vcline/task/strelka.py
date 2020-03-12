#!/usr/bin/env python

from itertools import product
from math import floor
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMs
from .base import ShellTask
from .manta import CallStructualVariantsWithManta
from .ref import CreateEvaluationIntervalListBED, FetchReferenceFASTA


@requires(PrepareCRAMs, FetchReferenceFASTA, CallStructualVariantsWithManta,
          CreateEvaluationIntervalListBED)
class CallVariantsWithStrelka(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['strelka_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0]]
                        )
                    ).joinpath(f'results/variants/somatic.{v}{i}')
                )
            )
            for v, i in product(['snvs.vcf.gz', 'indels.vcf.gz'], ['', '.tbi'])
        ]

    def run(self):
        output_file_paths = [o.path for o in self.output()]
        run_dir_path = str(Path(output_file_paths[0]).parent.parent.parent)
        run_id = Path(run_dir_path).name
        self.print_log(f'Call somatic variants with Strelka:\t{run_id}')
        config_script = self.cf['configureStrelkaSomaticWorkflow.py']
        run_script = str(Path(run_dir_path).joinpath('runWorkflow.py'))
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 1)
        input_cram_paths = [i[0].path for i in self.input()[0]]
        fa_path = self.input()[1].path
        manta_indel_vcf_path = str(
            Path(self.input()[2][0].path).joinpath(
                'candidateSmallIndels.vcf.gz'
            )
        )
        bed_path = self.input()[3][0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=config_script, cwd=self.cf['strelka_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {config_script}'
                + f' --tumorBam={input_cram_paths[0]}'
                + f' --normalBam={input_cram_paths[1]}'
                + f' --referenceFasta={fa_path}'
                + f' --indelCandidates={manta_indel_vcf_path}'
                + f' --callRegions={bed_path}'
                + f' --runDir={run_dir_path}'
            ),
            input_files_or_dirs=[
                *input_cram_paths, fa_path, manta_indel_vcf_path, bed_path
            ],
            output_files_or_dirs=[run_script, run_dir_path]
        )
        self.run_shell(
            args=(
                f'set -e && {run_script}'
                + f' --jobs={n_cpu}'
                + f' --memGb={memory_gb}'
                + ' --mode=local'
            ),
            input_files_or_dirs=[
                run_script, *input_cram_paths, fa_path, manta_indel_vcf_path,
                bed_path
            ],
            output_files_or_dirs=[*output_file_paths, run_dir_path]
        )


if __name__ == '__main__':
    luigi.run()
