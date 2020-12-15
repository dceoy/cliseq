#!/usr/bin/env python

from itertools import product
from math import floor
from pathlib import Path

import luigi
from ftarc.task.base import ShellTask
from ftarc.task.resource import FetchReferenceFASTA
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .ref import CreateEvaluationIntervalListBED


@requires(PrepareCRAMTumor, PrepareCRAMNormal, FetchReferenceFASTA,
          CreateEvaluationIntervalListBED)
class CallStructualVariantsWithManta(ShellTask):
    cf = luigi.DictParameter()
    priority = 40

    def output(self):
        run_dir = Path(self.cf['somatic_sv_manta_dir_path']).joinpath(
            create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(f'{run_dir.name}.manta.{v}SV.vcf.gz{s}')
            ) for v, s in product(['somatic', 'diploid'], ['', '.tbi'])
        ]

    def run(self):
        output_link_paths = [o.path for o in self.output()]
        run_dir = Path(output_link_paths[0]).parent
        run_id = run_dir.name
        self.print_log(f'Call somatic SVs with Manta:\t{run_id}')
        config_script = self.cf['configManta.py']
        run_script = run_dir.joinpath('runWorkflow.py')
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        python2 = self.cf['python2']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 4)
        input_cram_paths = [i[0].path for i in self.input()[0:2]]
        fa_path = self.input()[2][0].path
        bed_path = self.input()[3][0].path
        result_files = [
            run_dir.joinpath(f'results/variants/{v}.vcf.gz{s}')
            for v, s in product(
                [
                    'somaticSV', 'diploidSV', 'candidateSV',
                    'candidateSmallIndels'
                ],
                ['', '.tbi']
            )
        ]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[python2, config_script], cwd=run_dir,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {config_script}'
                + f' --tumorBam={input_cram_paths[0]}'
                + f' --normalBam={input_cram_paths[1]}'
                + f' --referenceFasta={fa_path}'
                + f' --callRegions={bed_path}'
                + f' --runDir={run_dir}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[*input_cram_paths, fa_path, bed_path],
            output_files_or_dirs=[run_script, run_dir]
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {run_script}'
                + f' --jobs={n_cpu}'
                + f' --memGb={memory_gb}'
                + ' --mode=local'
            ),
            input_files_or_dirs=[
                run_script, *input_cram_paths, fa_path, bed_path
            ],
            output_files_or_dirs=[*result_files, run_dir]
        )
        for p in output_link_paths:
            f = run_dir.joinpath('results/variants').joinpath(
                Path(p).name.split('.manta.')[-1]
            ).relative_to(run_dir)
            self.run_shell(args=f'ln -s {f} {p}', output_files_or_dirs=p)


if __name__ == '__main__':
    luigi.run()
