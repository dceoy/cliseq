#!/usr/bin/env python

from itertools import product
from math import floor
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .base import ShellTask
from .ref import CreateEvaluationIntervalListBED, FetchReferenceFASTA


@requires(PrepareCRAMTumor, PrepareCRAMNormal, FetchReferenceFASTA,
          CreateEvaluationIntervalListBED)
class CallStructualVariantsWithManta(ShellTask):
    cf = luigi.DictParameter()
    result_file_names = [
        (v + s) for v, s in product(
            [
                'somaticSV.vcf.gz', 'diploidSV.vcf.gz', 'candidateSV.vcf.gz',
                'candidateSmallIndels.vcf.gz'
            ],
            ['', '.tbi']
        )
    ]
    priority = 40

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['somatic_sv_manta_dir_path']).joinpath(
                    create_matched_id(*[i[0].path for i in self.input()[0:2]])
                    + s
                )
            ) for s in [
                '', '.manta.somaticSV.vcf.gz', '.manta.somaticSV.vcf.gz.tbi',
                '.manta.diploidSV.vcf.gz', '.manta.diploidSV.vcf.gz.tbi'
            ]
        ]

    def run(self):
        run_dir_path = self.output()[0].path
        run_id = Path(run_dir_path).name
        self.print_log(f'Call somatic SVs with Manta:\t{run_id}')
        config_script = self.cf['configManta.py']
        root_dir_path = self.cf['somatic_sv_manta_dir_path']
        run_script = str(Path(run_dir_path).joinpath('runWorkflow.py'))
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        python2 = self.cf['python2']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 1)
        input_cram_paths = [i[0].path for i in self.input()[0:2]]
        fa_path = self.input()[2][0].path
        bed_path = self.input()[3][0].path
        result_file_paths = [
            str(Path(run_dir_path).joinpath(f'results/variants/{n}'))
            for n in self.result_file_names
        ]
        output_link_paths = [o.path for o in self.output()[1:]]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[python2, config_script], cwd=root_dir_path,
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {config_script}'
                + f' --tumorBam={input_cram_paths[0]}'
                + f' --normalBam={input_cram_paths[1]}'
                + f' --referenceFasta={fa_path}'
                + f' --callRegions={bed_path}'
                + f' --runDir={run_dir_path}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[*input_cram_paths, fa_path, bed_path],
            output_files_or_dirs=[run_script, run_dir_path]
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
            output_files_or_dirs=[*result_file_paths, run_dir_path]
        )
        for l in output_link_paths:
            f = Path(run_dir_path).joinpath('results/variants').joinpath(
                Path(l).name.split('.manta.')[-1]
            ).relative_to(root_dir_path)
            self.run_shell(args=f'ln -s {f} {l}', output_files_or_dirs=l)


if __name__ == '__main__':
    luigi.run()
