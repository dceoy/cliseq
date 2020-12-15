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
from .bcftools import bcftools_concat_and_index
from .manta import CallStructualVariantsWithManta
from .ref import CreateEvaluationIntervalListBED


@requires(PrepareCRAMTumor, PrepareCRAMNormal, FetchReferenceFASTA,
          CreateEvaluationIntervalListBED, CallStructualVariantsWithManta)
class CallSomaticVariantsWithStrelka(ShellTask):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        run_dir = Path(self.cf['somatic_snv_indel_strelka_dir_path']).joinpath(
            create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(f'{run_dir.name}.strelka.somatic.vcf.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_dir = Path(output_vcf_path).parent
        run_id = run_dir.name
        self.print_log(f'Call somatic variants with Strelka:\t{run_id}')
        config_script = self.cf['configureStrelkaSomaticWorkflow.py']
        run_script = run_dir.joinpath('runWorkflow.py')
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        python2 = self.cf['python2']
        bcftools = self.cf['bcftools']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_mb = self.cf['memory_mb_per_worker']
        memory_gb = max(floor(memory_mb / 1024), 4)
        input_cram_paths = [i[0].path for i in self.input()[0:2]]
        fa_path = self.input()[2][0].path
        bed_path = self.input()[3][0].path
        manta_indel_vcf_path = Path(self.input()[4][0].path).parent.joinpath(
            'results/variants/candidateSmallIndels.vcf.gz'
        )
        result_files = [
            run_dir.joinpath(f'results/variants/somatic.{v}.vcf.gz{s}')
            for v, s in product(['snvs', 'indels'], ['', '.tbi'])
        ]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[python2, config_script, bcftools], cwd=run_dir,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {config_script}'
                + f' --tumorBam={input_cram_paths[0]}'
                + f' --normalBam={input_cram_paths[1]}'
                + f' --referenceFasta={fa_path}'
                + f' --indelCandidates={manta_indel_vcf_path}'
                + f' --callRegions={bed_path}'
                + f' --runDir={run_dir}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[
                *input_cram_paths, fa_path, manta_indel_vcf_path, bed_path
            ],
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
                run_script, *input_cram_paths, fa_path, manta_indel_vcf_path,
                bed_path
            ],
            output_files_or_dirs=[*result_files, run_dir]
        )
        bcftools_concat_and_index(
            shelltask=self, bcftools=bcftools,
            input_vcf_paths=[
                str(p) for p in result_files if p.name.endswith('.vcf.gz')
            ],
            output_vcf_path=output_vcf_path, n_cpu=n_cpu, memory_mb=memory_mb
        )


@requires(PrepareCRAMNormal, FetchReferenceFASTA,
          CreateEvaluationIntervalListBED)
class CallGermlineVariantsWithStrelka(ShellTask):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        run_dir = Path(
            self.cf['germline_snv_indel_strelka_dir_path']
        ).joinpath(Path(self.input()[0][0].path).stem)
        return [
            luigi.LocalTarget(
                run_dir.joinpath(
                    f'{run_dir.name}.strelka.germline.vcf.gz{s}'
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_link_paths = [o.path for o in self.output()]
        run_dir = Path(output_link_paths[0]).parent
        run_id = run_dir.name
        self.print_log(f'Call germline variants with Strelka:\t{run_id}')
        config_script = self.cf['configureStrelkaGermlineWorkflow.py']
        run_script = run_dir.joinpath('runWorkflow.py')
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 4)
        input_cram_path = self.input()[0][0].path
        fa_path = self.input()[1][0].path
        bed_path = self.input()[2][0].path
        result_files = [
            run_dir.joinpath(f'results/variants/{v}.vcf.gz{s}')
            for v, s in product(['variants', 'genome'], ['', '.tbi'])
        ]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=config_script, cwd=run_dir,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {config_script}'
                + f' --bam={input_cram_path}'
                + f' --referenceFasta={fa_path}'
                + f' --callRegions={bed_path}'
                + f' --runDir={run_dir}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[input_cram_path, fa_path, bed_path],
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
                run_script, input_cram_path, fa_path, bed_path
            ],
            output_files_or_dirs=[*result_files, run_dir]
        )
        for p in output_link_paths:
            f = run_dir.joinpath('results/variants').joinpath(
                'variants.' + Path(p).name.split('.strelka.germline.')[-1]
            ).relative_to(run_dir)
            self.run_shell(args=f'ln -s {f} {p}', output_files_or_dirs=p)


if __name__ == '__main__':
    luigi.run()
