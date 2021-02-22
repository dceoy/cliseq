#!/usr/bin/env python

import os
from itertools import product
from math import floor
from pathlib import Path

import luigi
from ftarc.task.resource import FetchReferenceFasta
from luigi.util import requires

from .core import VclineTask
from .cram import PrepareCramNormal, PrepareCramTumor
from .manta import CallStructualVariantsWithManta
from .resource import CreateEvaluationIntervalListBed


@requires(PrepareCramTumor, PrepareCramNormal, FetchReferenceFasta,
          CreateEvaluationIntervalListBed, CallStructualVariantsWithManta)
class CallSomaticVariantsWithStrelka(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.cf['somatic_snv_indel_strelka_dir_path']).joinpath(
            self.create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(f'{run_dir.name}.strelka.somatic.vcf.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf = Path(self.output()[0].path)
        run_dir = output_vcf.parent
        run_id = run_dir.name
        self.print_log(f'Call somatic variants with Strelka:\t{run_id}')
        config_script = Path(self.cf['configureStrelkaSomaticWorkflow.py'])
        run_script = run_dir.joinpath('runWorkflow.py')
        python2 = self.cf['python2']
        bcftools = self.cf['bcftools']
        pythonpath = '{0}:{1}'.format(
            Path(config_script).parent.parent.joinpath('lib/python'),
            (os.getenv('PYTHONPATH') or '')
        )
        memory_gb = max(floor(self.memory_mb / 1024), 4)
        input_crams = [Path(i[0].path) for i in self.input()[0:2]]
        fa = Path(self.input()[2][0].path)
        bed = Path(self.input()[3][0].path)
        manta_indel_vcf = Path(self.input()[4][0].path).parent.joinpath(
            'results/variants/candidateSmallIndels.vcf.gz'
        )
        result_files = [
            run_dir.joinpath(f'results/variants/somatic.{v}.vcf.gz{s}')
            for v, s in product(['snvs', 'indels'], ['', '.tbi'])
        ]
        self.setup_shell(
            run_id=run_id, commands=[python2, config_script, bcftools],
            cwd=run_dir, **self.sh_config, env={'PYTHONPATH': pythonpath}
        )
        self.run_shell(
            args=(
                f'set -e && {python2} {config_script}'
                + f' --tumorBam={input_crams[0]}'
                + f' --normalBam={input_crams[1]}'
                + f' --referenceFasta={fa}'
                + f' --indelCandidates={manta_indel_vcf}'
                + f' --callRegions={bed}'
                + f' --runDir={run_dir}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[*input_crams, fa, manta_indel_vcf, bed],
            output_files_or_dirs=[run_script, run_dir]
        )
        self.run_shell(
            args=(
                f'set -e && {python2} {run_script} --mode=local'
                + f' --jobs={self.n_cpu} --memGb={memory_gb}'
            ),
            input_files_or_dirs=[
                run_script, *input_crams, fa, manta_indel_vcf, bed
            ],
            output_files_or_dirs=[*result_files, run_dir]
        )
        self.bcftools_concat(
            input_vcf_paths=[
                o for o in result_files if o.name.endswith('.vcf.gz')
            ],
            output_vcf_path=output_vcf, bcftools=bcftools, n_cpu=self.n_cpu,
            memory_mb=self.memory_mb, index_vcf=True, remove_input=False
        )


@requires(PrepareCramNormal, FetchReferenceFasta,
          CreateEvaluationIntervalListBed)
class CallGermlineVariantsWithStrelka(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
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
        output_links = [Path(o.path) for o in self.output()]
        run_dir = output_links[0].parent
        run_id = run_dir.name
        self.print_log(f'Call germline variants with Strelka:\t{run_id}')
        config_script = Path(self.cf['configureStrelkaGermlineWorkflow.py'])
        run_script = run_dir.joinpath('runWorkflow.py')
        python2 = self.cf['python2']
        pythonpath = '{0}:{1}'.format(
            Path(config_script).parent.parent.joinpath('lib/python'),
            (os.getenv('PYTHONPATH') or '')
        )
        memory_gb = max(floor(self.memory_mb / 1024), 4)
        input_cram = Path(self.input()[0][0].path)
        fa = Path(self.input()[1][0].path)
        bed = Path(self.input()[2][0].path)
        result_files = [
            run_dir.joinpath(f'results/variants/{v}.vcf.gz{s}')
            for v, s in product(['variants', 'genome'], ['', '.tbi'])
        ]
        self.setup_shell(
            run_id=run_id, commands=[python2, config_script], cwd=run_dir,
            **self.sh_config, env={'PYTHONPATH': pythonpath}
        )
        self.run_shell(
            args=(
                f'set -e && {python2} {config_script}'
                + f' --bam={input_cram}'
                + f' --referenceFasta={fa}'
                + f' --callRegions={bed}'
                + f' --runDir={run_dir}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[input_cram, fa, bed],
            output_files_or_dirs=[run_script, run_dir]
        )
        self.run_shell(
            args=(
                f'set -e && {python2} {run_script} --mode=local'
                + f' --jobs={self.n_cpu} --memGb={memory_gb}'
            ),
            input_files_or_dirs=[run_script, input_cram, fa, bed],
            output_files_or_dirs=[*result_files, run_dir]
        )
        for o in output_links:
            f = run_dir.joinpath('results/variants').joinpath(
                'variants.' + o.name.split('.strelka.germline.')[-1]
            ).relative_to(run_dir)
            self.run_shell(args=f'ln -s {f} {o}', output_files_or_dirs=o)


if __name__ == '__main__':
    luigi.run()
