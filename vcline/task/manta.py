#!/usr/bin/env python

from itertools import product
from math import floor
from pathlib import Path

import luigi
from ftarc.task.resource import FetchReferenceFasta
from luigi.util import requires

from .core import VclineTask
from .cram import PrepareCramNormal, PrepareCramTumor
from .resource import CreateEvaluationIntervalListBed


@requires(PrepareCramTumor, PrepareCramNormal, FetchReferenceFasta,
          CreateEvaluationIntervalListBed)
class CallStructualVariantsWithManta(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 40

    def output(self):
        run_dir = Path(self.cf['somatic_sv_manta_dir_path']).joinpath(
            self.create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(f'{run_dir.name}.manta.{v}SV.vcf.gz{s}')
            ) for v, s in product(['somatic', 'diploid'], ['', '.tbi'])
        ]

    def run(self):
        output_links = [Path(o.path) for o in self.output()]
        run_dir = output_links[0].parent
        run_id = run_dir.name
        self.print_log(f'Call somatic SVs with Manta:\t{run_id}')
        config_script = Path(self.cf['configManta.py']).resolve()
        run_script = run_dir.joinpath('runWorkflow.py')
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        python2 = self.cf['python2']
        memory_gb = max(floor(self.memory_mb / 1024), 4)
        input_crams = [Path(i[0].path) for i in self.input()[0:2]]
        fa = Path(self.input()[2][0].path)
        bed = Path(self.input()[3][0].path)
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
            run_id=run_id, commands=[python2, config_script], cwd=run_dir,
            **self.sh_config, env={'PYTHONPATH': str(pythonpath)}
        )
        self.run_shell(
            args=(
                f'set -e && {config_script}'
                + f' --tumorBam={input_crams[0]}'
                + f' --normalBam={input_crams[1]}'
                + f' --referenceFasta={fa}'
                + f' --callRegions={bed}'
                + f' --runDir={run_dir}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[*input_crams, fa, bed],
            output_files_or_dirs=[run_script, run_dir]
        )
        self.run_shell(
            args=(
                f'set -e && {run_script} --mode=local'
                + f' --jobs={self.n_cpu} --memGb={memory_gb}'
            ),
            input_files_or_dirs=[run_script, *input_crams, fa, bed],
            output_files_or_dirs=[*result_files, run_dir]
        )
        for o in output_links:
            f = run_dir.joinpath('results/variants').joinpath(
                o.name.split('.manta.')[-1]
            ).relative_to(run_dir)
            self.run_shell(args=f'ln -s {f} {o}', output_files_or_dirs=o)


if __name__ == '__main__':
    luigi.run()
