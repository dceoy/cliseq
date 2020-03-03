#!/usr/bin/env python

from itertools import product
from math import floor
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMs
from .base import ShellTask
from .ref import CreateEvaluationIntervalListBED, FetchReferenceFASTA


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateEvaluationIntervalListBED)
class CallStructualVariantsWithManta(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['manta_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0]]
                        )
                    ).joinpath(f'results/variants/{v}{i}')
                )
            ) for v, i in product(
                [
                    'somaticSV.vcf.gz', 'diploidSV.vcf.gz', 'tumorSV.vcf.gz',
                    'candidateSV.vcf.gz', 'candidateSmallIndels.vcf.gz'
                ],
                ['', '.tbi']
            )
        ]

    def run(self):
        output_file_paths = [o.path for o in self.output()]
        run_dir_path = str(Path(output_file_paths[0]).parent.parent.parent)
        run_id = Path(run_dir_path).name
        self.print_log(f'Call somatic SVs with Manta:\t{run_id}')
        config_script = self.cf['configManta.py']
        run_script = str(Path(run_dir_path).joinpath('runWorkflow.py'))
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 1)
        input_cram_paths = [i[0].path for i in self.input()[0]]
        fa_path = self.input()[1].path
        bed_path = self.input()[2][0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=config_script, cwd=self.cf['manta_dir_path']
        )
        self.run_shell(
            args=(
                f'set -e && {config_script}'
                + f' --tumorBam={input_cram_paths[0]}'
                + f' --normalBam={input_cram_paths[1]}'
                + f' --referenceFasta={fa_path}'
                + f' --callRegions={bed_path}'
                + f' --runDir={run_dir_path}'
            ),
            input_files=[*input_cram_paths, fa_path, bed_path],
            output_files=run_script
        )
        self.run_shell(
            args=(
                f'set -e && {run_script}'
                + f' --jobs={n_cpu}'
                + f' --memGb={memory_gb}'
                + ' --mode=local'
            ),
            input_files=[
                run_script, *input_cram_paths, fa_path, bed_path
            ],
            output_files=output_file_paths
        )


if __name__ == '__main__':
    luigi.run()
