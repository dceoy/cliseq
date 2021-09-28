#!/usr/bin/env python

import re
import sys
from itertools import product
from pathlib import Path

import luigi
from luigi.util import requires

from .core import VclineTask
from .cram import PrepareCramNormal, PrepareCramTumor
from .resource import CreateEvaluationIntervalListBed, FetchReferenceFasta


@requires(CreateEvaluationIntervalListBed, FetchReferenceFasta)
class CreateExclusionIntervalListBed(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        return [
            luigi.LocalTarget(p + s) for p, s in product(
                [
                    re.sub(
                        r'\.bed\.gz$', '.exclusion.bed.gz',
                        self.input()[0][0].path
                    ),
                    (self.input()[1][0].path + '.bed.gz')
                ],
                ['', '.tbi']
            )
        ]

    def run(self):
        evaluation_bed = Path(self.input()[0][0].path)
        run_id = Path(evaluation_bed.stem).stem
        self.print_log(f'Create an exclusion interval_list BED:\t{run_id}')
        bedtools = self.cf['bedtools']
        bgzip = self.cf['bgzip']
        tabix = self.cf['tabix']
        fai_path = self.input()[1][1].path
        exclusion_bed_path = self.output()[0].path
        genome_bed_path = self.output()[2].path
        self.setup_shell(
            run_id=run_id, commands=[bgzip, tabix], cwd=evaluation_bed.parent,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {sys.executable}'
                + ' -c \'{}\''.format(
                    'from fileinput import input; '
                    '[print("{0}\\t0\\t{1}".format(*s.split()[:2]))'
                    ' for s in input()];'
                ) + f' {fai_path}'
                + f' | {bgzip} -@ {self.n_cpu} -c > {genome_bed_path}'
            ),
            input_files_or_dirs=fai_path, output_files_or_dirs=genome_bed_path
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {bedtools} subtract'
                + f' -a {genome_bed_path} -b {evaluation_bed}'
                + f' | {bgzip} -@ {self.n_cpu} -c > {exclusion_bed_path}'
            ),
            input_files_or_dirs=[genome_bed_path, evaluation_bed],
            output_files_or_dirs=exclusion_bed_path
        )
        for p in [genome_bed_path, exclusion_bed_path]:
            self.tabix_tbi(tsv_path=p, tabix=tabix, preset='bed')


@requires(PrepareCramTumor, PrepareCramNormal, FetchReferenceFasta,
          CreateExclusionIntervalListBed)
class CallSomaticStructualVariantsWithDelly(VclineTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.cf['somatic_sv_delly_dir_path']).joinpath(
            self.create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(run_dir.name + f'.delly.{s}')
            ) for s in [
                'filtered.vcf.gz', 'filtered.vcf.gz.tbi', 'bcf', 'bcf.csi',
                'filtered.bcf', 'filtered.bcf.csi'
            ]
        ]

    def run(self):
        output_vcf = Path(self.output()[0].path)
        run_dir = output_vcf.parent
        run_id = run_dir.name
        self.print_log(f'Call somatic SVs with Delly:\t{run_id}')
        delly = self.cf['delly']
        bcftools = self.cf['bcftools']
        input_crams = [Path(i[0].path) for i in self.input()[0:2]]
        fa = Path(self.input()[2][0].path)
        exclusion_bed = Path(self.input()[3][0].path)
        raw_bcf = Path(self.output()[2].path)
        filtered_bcf = Path(self.output()[4].path)
        samples_tsv = run_dir.joinpath(f'{raw_bcf.stem}.tsv')
        self.setup_shell(
            run_id=run_id, commands=[delly, bcftools], cwd=run_dir,
            **self.sh_config, env={'OMP_NUM_THREADS': str(self.n_cpu)}
        )
        self.run_shell(
            args=(
                f'set -e && {delly} call'
                + f' --genome {fa}'
                + f' --exclude {exclusion_bed}'
                + f' --outfile {raw_bcf}'
                + ''.join(f' {p}' for p in input_crams)
            ),
            input_files_or_dirs=[*input_crams, fa, exclusion_bed],
            output_files_or_dirs=[raw_bcf, f'{raw_bcf}.csi', run_dir]
        )
        self.run_shell(
            args=(
                'set -e && echo -e'
                + ' "{0}\\ttumor\\n{1}\\tnormal\\n"'.format(*self.sample_names)
                + f' | tee {samples_tsv}'
            ),
            output_files_or_dirs=samples_tsv
        )
        self.run_shell(
            args=(
                f'set -e && {delly} filter --filter somatic'
                + f' --samples {samples_tsv}'
                + f' --outfile {filtered_bcf}'
                + f' {raw_bcf}'
            ),
            input_files_or_dirs=[raw_bcf, samples_tsv],
            output_files_or_dirs=[filtered_bcf, f'{filtered_bcf}.csi']
        )
        self.bcftools_sort(
            input_vcf_path=filtered_bcf, output_vcf_path=output_vcf,
            bcftools=bcftools, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            index_vcf=True, remove_input=False
        )


@requires(PrepareCramNormal, FetchReferenceFasta,
          CreateExclusionIntervalListBed)
class CallGermlineStructualVariantsWithDelly(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.cf['germline_sv_delly_dir_path']).joinpath(
            Path(self.input()[0][0].path).stem
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(run_dir.name + f'.delly.{s}')
            ) for s in ['bcf', 'bcf.csi', 'vcf.gz', 'vcf.gz.tbi']
        ]

    def run(self):
        output_bcf = Path(self.output()[0].path)
        run_dir = output_bcf.parent
        run_id = run_dir.name
        self.print_log(f'Call germline SVs with Delly:\t{run_id}')
        delly = self.cf['delly']
        bcftools = self.cf['bcftools']
        input_cram = Path(self.input()[0][0].path)
        fa = Path(self.input()[2][0].path)
        exclusion_bed = Path(self.input()[3][0].path)
        output_vcf = Path(self.output()[2].path)
        self.setup_shell(
            run_id=run_id, commands=[delly, bcftools], cwd=run_dir,
            **self.sh_config, env={'OMP_NUM_THREADS': str(self.n_cpu)}
        )
        self.run_shell(
            args=(
                f'set -e && {delly} call'
                + f' --genome {fa}'
                + f' --exclude {exclusion_bed}'
                + f' --outfile {output_bcf}'
                + f' {input_cram}'
            ),
            input_files_or_dirs=[input_cram, fa, exclusion_bed],
            output_files_or_dirs=[output_bcf, f'{output_bcf}.csi', run_dir]
        )
        self.bcftools_sort(
            input_vcf_path=output_bcf, output_vcf_path=output_vcf,
            bcftools=bcftools, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            index_vcf=True, remove_input=False
        )


if __name__ == '__main__':
    luigi.run()
