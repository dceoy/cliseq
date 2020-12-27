#!/usr/bin/env python

import re
from itertools import chain
from pathlib import Path

import luigi
from ftarc.task.picard import CreateSequenceDictionary
from ftarc.task.resource import FetchReferenceFasta
from luigi.util import requires

from .core import VclineTask
from .cram import PrepareCramNormal, PrepareCramTumor
from .haplotypecaller import SplitEvaluationIntervals
from .resource import (CreateGnomadBiallelicSnpVcf,
                       FetchEvaluationIntervalList, FetchGnomadVcf)


class GetPileupSummaries(VclineTask):
    cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    gnomad_common_biallelic_vcf_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.cram_path).stem + '.pileup.table'
            )
        )

    def run(self):
        cram = Path(self.cram_path).resolve()
        run_id = cram.stem
        self.print_log(f'Get pileup summary:\t{run_id}')
        output_pileup_table = Path(self.output().path)
        fa = Path(self.fa_path).resolve()
        evaluation_interval = Path(self.evaluation_interval_path).resolve()
        gnomad_common_biallelic_vcf = Path(
            self.gnomad_common_biallelic_vcf_path
        ).resolve()
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=output_pileup_table.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} GetPileupSummaries'
                + f' --input {cram}'
                + f' --reference {fa}'
                + f' --variant {gnomad_common_biallelic_vcf}'
                + f' --intervals {evaluation_interval}'
                + f' --output {output_pileup_table}'
                + ' --disable-bam-index-caching '
                + str(self.save_memory).lower()
            ),
            input_files_or_dirs=[
                cram, fa, evaluation_interval, gnomad_common_biallelic_vcf
            ],
            output_files_or_dirs=output_pileup_table
        )


@requires(PrepareCramTumor, PrepareCramNormal, FetchReferenceFasta,
          FetchEvaluationIntervalList, CreateGnomadBiallelicSnpVcf,
          CreateSequenceDictionary)
class CalculateContamination(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        output_path_prefix = Path(self.cf['qc_dir_path']).joinpath(
            'contamination'
        ).joinpath(
            self.create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(f'{output_path_prefix}.{s}.table')
            for s in ['contamination', 'segment']
        ]

    def run(self):
        output_contamination_table = Path(self.output()[0].path)
        run_dir = output_contamination_table.parent
        gatk = self.cf['gatk']
        input_targets = yield [
            GetPileupSummaries(
                cram_path=self.input()[i][0].path,
                fa_path=self.input()[2][0].path,
                evaluation_interval_path=self.input()[3].path,
                gnomad_common_biallelic_vcf_path=self.input()[4][0].path,
                dest_dir_path=str(run_dir), gatk=gatk,
                save_memory=self.cf['save_memory']
            ) for i in range(2)
        ]
        run_id = '.'.join(output_contamination_table.name.split('.')[:-2])
        self.print_log(f'Calculate cross-sample contamination:\t{run_id}')
        pileup_tables = [Path(i.path) for i in input_targets]
        output_segment_table = Path(self.output()[1].path)
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=run_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} CalculateContamination'
                + f' --input {pileup_tables[0]}'
                + f' --matched-normal {pileup_tables[1]}'
                + f' --output {output_contamination_table}'
                + f' --tumor-segmentation {output_segment_table}'
            ),
            input_files_or_dirs=pileup_tables,
            output_files_or_dirs=[
                output_contamination_table, output_segment_table
            ]
        )


@requires(PrepareCramTumor, PrepareCramNormal, FetchReferenceFasta,
          SplitEvaluationIntervals, FetchGnomadVcf,
          CreateSequenceDictionary)
class CallVariantsWithMutect2(VclineTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        run_dir = Path(self.cf['somatic_snv_indel_gatk_dir_path']).joinpath(
            self.create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(run_dir.joinpath(f'{run_dir.name}.mutect2.{s}'))
            for s in [
                'vcf.gz', 'vcf.gz.tbi', 'vcf.gz.stats', 'cram', 'cram.crai',
                'read-orientation-model.tar.gz'
            ]
        ]

    def run(self):
        output_vcf = Path(self.output()[0].path)
        intervals = [Path(i.path) for i in self.input()[3]]
        skip_interval_split = (len(intervals) == 1)
        fa = Path(self.input()[2][0].path)
        input_crams = [Path(i[0].path) for i in self.input()[0:2]]
        gnomad_vcf = Path(self.input()[4][0].path)
        output_path_prefix = '.'.join(str(output_vcf).split('.')[:-2])
        if skip_interval_split:
            tmp_prefixes = [output_path_prefix]
        else:
            tmp_prefixes = [
                '{0}.{1}'.format(output_path_prefix, o.stem) for o in intervals
            ]
        input_targets = yield [
            Mutect2(
                input_cram_paths=input_crams, fa_path=fa,
                gnomad_vcf_path=gnomad_vcf, evaluation_interval_path=str(o),
                normal_name=self.sample_names[1], output_path_prefix=s,
                cf=self.cf
            ) for o, s in zip(intervals, tmp_prefixes)
        ]
        run_id = '.'.join(output_vcf.name.split('.')[:-3])
        self.print_log(f'Call somatic variants with Mutect2:\t{run_id}')
        output_stats = Path(self.output()[2].path)
        output_cram = Path(self.output()[3].path)
        ob_priors = Path(self.output()[5].path)
        f1r2s = [f'{s}.f1r2.tar.gz' for s in tmp_prefixes]
        gatk = self.cf['gatk']
        samtools = self.cf['samtools']
        self.setup_shell(
            run_id=run_id, commands=[gatk, samtools], cwd=output_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} LearnReadOrientationModel'
                + ''.join([f' --input {f}' for f in f1r2s])
                + f' --output {ob_priors}'
            ),
            input_files_or_dirs=f1r2s, output_files_or_dirs=ob_priors
        )
        if skip_interval_split:
            tmp_bam = f'{tmp_prefixes[0]}.bam'
            self.samtools_view(
                input_sam_path=tmp_bam, fa_path=fa,
                output_sam_path=output_cram, samtools=samtools,
                n_cpu=self.n_cpu, index_sam=True, remove_input=True
            )
        else:
            tmp_vcfs = [Path(f'{s}.vcf.gz') for s in tmp_prefixes]
            self.picard_mergevcfs(
                input_vcf_paths=tmp_vcfs, output_vcf_path=output_vcf,
                picard=gatk, remove_input=False
            )
            self.run_shell(
                args=(
                    f'set -e && {gatk} MergeVcfs'
                    + ''.join([f' --INPUT {v}' for v in tmp_vcfs])
                    + f' --OUTPUT {output_vcf}'
                ),
                input_files_or_dirs=tmp_vcfs,
                output_files_or_dirs=[output_vcf, f'{output_vcf}.tbi']
            )
            tmp_statses = [Path(f'{s}.vcf.gz.stats') for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {gatk} MergeMutectStats'
                    + ''.join([f' --stats {s}' for s in tmp_statses])
                    + f' --output {output_stats}'
                ),
                input_files_or_dirs=tmp_statses,
                output_files_or_dirs=output_stats
            )
            self.samtools_merge(
                input_sam_paths=[f'{s}.bam' for s in tmp_prefixes],
                fa_path=fa, output_sam_path=output_cram, samtools=samtools,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb, index_sam=True,
                remove_input=False
            )
            self.remove_files_and_dirs(
                *chain.from_iterable(
                    [o.path for o in t] for t in input_targets
                )
            )


class Mutect2(VclineTask):
    input_cram_paths = luigi.ListParameter()
    fa_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    normal_name = luigi.Parameter()
    output_path_prefix = luigi.Parameter()
    cf = luigi.DictParameter()
    message = luigi.Parameter(default='')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(f'{self.output_path_prefix}.{s}') for s in [
                'vcf.gz', 'vcf.gz.tbi', 'vcf.gz.stats', 'bam', 'f1r2.tar.gz'
            ]
        ]

    def run(self):
        if self.message:
            self.print_log(self.message)
        input_crams = [Path(p).resolve() for p in self.input_cram_paths]
        fa = Path(self.fa_path).resolve()
        gnomad_vcf = Path(self.gnomad_vcf_path).resolve()
        evaluation_interval = Path(self.evaluation_interval_path).resolve()
        output_files = [Path(o.path) for o in self.output()]
        output_vcf = output_files[0]
        run_dir = output_vcf.parent
        gatk = self.cf['gatk']
        self.setup_shell(
            run_id='.'.join(output_vcf.name.split('.')[:-2]), commands=gatk,
            cwd=run_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} Mutect2'
                + ''.join([f' --input {c}' for c in input_crams])
                + f' --reference {fa}'
                + f' --intervals {evaluation_interval}'
                + f' --germline-resource {gnomad_vcf}'
                + f' --output {output_vcf}'
                + f' --bam-output {output_files[3]}'
                + f' --f1r2-tar-gz {output_files[4]}'
                + f' --normal-sample {self.normal_name}'
                + ' --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP'
                + f' --native-pair-hmm-threads {self.n_cpu}'
                + ' --max-mnp-distance 0'
                + ' --create-output-bam-index false'
                + ' --disable-bam-index-caching '
                + str(self.cf['save_memory']).lower()
            ),
            input_files_or_dirs=[
                *input_crams, fa, evaluation_interval, gnomad_vcf
            ],
            output_files_or_dirs=[*output_files, run_dir]
        )


@requires(CallVariantsWithMutect2, FetchReferenceFasta,
          FetchEvaluationIntervalList, CalculateContamination,
          CreateSequenceDictionary)
class FilterMutectCalls(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        output_vcf_path = re.sub(
            r'\.vcf\.gz$', '.filtered.vcf.gz', self.input()[0][0].path
        )
        return [
            luigi.LocalTarget(output_vcf_path + s)
            for s in ['', '.tbi', '.stats']
        ]

    def run(self):
        output_filtered_vcf = Path(self.output()[0].path)
        run_id = '.'.join(output_filtered_vcf.name.split('.')[:-4])
        self.print_log(f'Filter somatic variants called by Mutect2:\t{run_id}')
        mutect_vcf = Path(self.input()[0][0].path)
        mutect_stats = Path(self.input()[0][2].path)
        ob_priors = Path(self.input()[0][5].path)
        fa = Path(self.input()[1][0].path)
        evaluation_interval = Path(self.input()[2].path)
        contamination_table = Path(self.input()[3][0].path)
        segment_table = Path(self.input()[3][1].path)
        output_filtering_stats = Path(self.output()[2].path)
        gatk = self.cf['gatk']
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=output_filtered_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} FilterMutectCalls'
                + f' --reference {fa}'
                + f' --intervals {evaluation_interval}'
                + f' --variant {mutect_vcf}'
                + f' --stats {mutect_stats}'
                + f' --contamination-table {contamination_table}'
                + f' --tumor-segmentation {segment_table}'
                + f' --orientation-bias-artifact-prior {ob_priors}'
                + f' --output {output_filtered_vcf}'
                + f' --filtering-stats {output_filtering_stats}'
            ),
            input_files_or_dirs=[
                mutect_vcf, fa, evaluation_interval, mutect_stats, ob_priors,
                contamination_table, segment_table,
            ],
            output_files_or_dirs=[output_filtered_vcf, output_filtering_stats]
        )


if __name__ == '__main__':
    luigi.run()
