#!/usr/bin/env python

import re
from itertools import chain
from pathlib import Path

import luigi
from ftarc.task.picard import CreateSequenceDictionary
from ftarc.task.resource import (FetchDbsnpVcf, FetchMillsIndelVcf,
                                 FetchReferenceFasta)
from luigi.util import requires

from .core import VclineTask
from .cram import PrepareCramNormal
from .resource import FetchEvaluationIntervalList, FetchHapmapVcf


@requires(FetchEvaluationIntervalList, FetchReferenceFasta,
          CreateSequenceDictionary)
class SplitEvaluationIntervals(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        if self.cf['n_worker'] > 1:
            run_dir = Path(self.cf['qc_dir_path']).joinpath(
                'intervals/{0}.split_in_{1}'.format(
                    Path(self.input()[0].path).stem, self.cf['n_worker']
                )
            )
            return [
                luigi.LocalTarget(
                    run_dir.joinpath(f'{i:04d}-scattered.interval_list')
                ) for i in range(self.cf['n_worker'])
            ]
        else:
            return [luigi.LocalTarget(self.input()[0].path)]

    def run(self):
        input_interval = Path(self.input()[0].path)
        run_id = input_interval.stem
        output_intervals = [Path(o.path) for o in self.output()]
        scatter_count = len(output_intervals)
        self.print_log(f'Split an interval list in {scatter_count}:\t{run_id}')
        fa = Path(self.input()[1][0].path)
        run_dir = output_intervals[0].parent
        gatk = self.cf['gatk']
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
                f'set -e && {gatk} SplitIntervals'
                + f' --reference {fa}'
                + f' --intervals {input_interval}'
                + f' --scatter-count {scatter_count}'
                + f' --output {run_dir}'
            ),
            input_files_or_dirs=[input_interval, fa],
            output_files_or_dirs=[*output_intervals, run_dir]
        )


@requires(PrepareCramNormal, FetchReferenceFasta, FetchDbsnpVcf,
          SplitEvaluationIntervals, CreateSequenceDictionary)
class CallVariantsWithHaplotypeCaller(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        run_dir = Path(self.cf['germline_snv_indel_gatk_dir_path']).joinpath(
            Path(self.input()[0][0].path).stem
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(f'{run_dir.name}.haplotypecaller.{s}')
            ) for s in ['g.vcf.gz', 'g.vcf.gz.tbi', 'cram', 'cram.crai']
        ]

    def run(self):
        output_gvcf = Path(self.output()[0].path)
        intervals = [Path(i.path) for i in self.input()[3]]
        skip_interval_split = (len(intervals) == 1)
        fa = Path(self.input()[1][0].path)
        input_cram = Path(self.input()[0][0].path)
        dbsnp_vcf = Path(self.input()[2][0].path)
        output_path_prefix = '.'.join(str(output_gvcf).split('.')[:-3])
        if skip_interval_split:
            tmp_prefixes = [output_path_prefix]
        else:
            tmp_prefixes = [
                '{0}.{1}'.format(output_path_prefix, o.stem) for o in intervals
            ]
        input_targets = yield [
            HaplotypeCaller(
                input_cram_path=str(input_cram), fa_path=str(fa),
                dbsnp_vcf_path=str(dbsnp_vcf), evaluation_interval_path=str(o),
                output_path_prefix=s, cf=self.cf
            ) for o, s in zip(intervals, tmp_prefixes)
        ]
        run_id = '.'.join(output_gvcf.name.split('.')[:-4])
        self.print_log(
            f'Call germline variants with HaplotypeCaller:\t{run_id}'
        )
        output_cram = Path(self.output()[2].path)
        gatk = self.cf['gatk']
        samtools = self.cf['samtools']
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=output_gvcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        if skip_interval_split:
            tmp_bam = Path(f'{tmp_prefixes[0]}.bam')
            self.samtools_view(
                input_sam_path=tmp_bam, fa_path=fa,
                output_sam_path=output_cram, samtools=samtools,
                n_cpu=self.n_cpu, index_sam=True, remove_input=True
            )
        else:
            tmp_gvcfs = [Path(f'{s}.g.vcf.gz') for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {gatk} CombineGvcfs'
                    + f' --reference {fa}'
                    + ''.join([f' --variant {p}' for p in tmp_gvcfs])
                    + f' --output {output_gvcf}'
                ),
                input_files_or_dirs=[*tmp_gvcfs, fa],
                output_files_or_dirs=[output_gvcf, f'{output_gvcf}.tbi']
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


class HaplotypeCaller(VclineTask):
    input_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dbsnp_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    output_path_prefix = luigi.Parameter()
    cf = luigi.DictParameter()
    message = luigi.Parameter(default='')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(f'{self.output_path_prefix}.{s}')
            for s in ['g.vcf.gz', 'g.vcf.gz.tbi', 'bam']
        ]

    def run(self):
        if self.message:
            self.print_log(self.message)
        input_cram = Path(self.input_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        dbsnp_vcf = Path(self.dbsnp_vcf_path).resolve()
        evaluation_interval = Path(self.evaluation_interval_path).resolve()
        output_files = [Path(o.path) for o in self.output()]
        output_gvcf = output_files[0]
        run_dir = output_gvcf.parent
        gatk = self.cf['gatk']
        self.setup_shell(
            run_id='.'.join(output_gvcf.name.split('.')[:-3]), commands=gatk,
            cwd=run_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} HaplotypeCaller'
                + f' --input {input_cram}'
                + f' --reference {fa}'
                + f' --dbsnp {dbsnp_vcf}'
                + f' --intervals {evaluation_interval}'
                + f' --output {output_gvcf}'
                + f' --bam-output {output_files[2]}'
                + ' --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP'
                + f' --native-pair-hmm-threads {self.n_cpu}'
                + ' --emit-ref-confidence GVCF'
                + ''.join(
                    [
                        f' --annotation-group {g}' for g in [
                            'StandardAnnotation', 'AS_StandardAnnotation',
                            'StandardHCAnnotation'
                        ]
                    ] + [f' --gvcf-gq-bands {i}' for i in range(10, 100, 10)]
                )
                + ' --create-output-bam-index false'
                + ' --disable-bam-index-caching '
                + str(self.cf['save_memory']).lower()
            ),
            input_files_or_dirs=[
                input_cram, fa, dbsnp_vcf, evaluation_interval
            ],
            output_files_or_dirs=[*output_files, run_dir]
        )


@requires(CallVariantsWithHaplotypeCaller, FetchReferenceFasta,
          FetchDbsnpVcf, FetchEvaluationIntervalList, CreateSequenceDictionary)
class GenotypeHaplotypeCallerGvcf(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        output_vcf_path = re.sub(
            r'\.g\.vcf\.gz$', '.vcf.gz', self.input()[0][0].path
        )
        return [luigi.LocalTarget(output_vcf_path + s) for s in ['', '.tbi']]

    def run(self):
        output_vcf = Path(self.output()[0].path)
        run_id = '.'.join(output_vcf.name.split('.')[:-3])
        self.print_log(f'Genotype a HaplotypeCaller GVCF:\t{run_id}')
        gvcf = Path(self.input()[0][0].path)
        fa = Path(self.input()[1][0].path)
        dbsnp_vcf = Path(self.input()[2][0].path)
        evaluation_interval = Path(self.input()[3].path)
        gatk = self.cf['gatk']
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=output_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} GenotypeGvcfs'
                + f' --reference {fa}'
                + f' --variant {gvcf}'
                + f' --dbsnp {dbsnp_vcf}'
                + f' --intervals {evaluation_interval}'
                + f' --output {output_vcf}'
                + ' --disable-bam-index-caching '
                + str(self.cf['save_memory']).lower()
            ),
            input_files_or_dirs=[gvcf, fa, dbsnp_vcf, evaluation_interval],
            output_files_or_dirs=[output_vcf, f'{output_vcf}.tbi']
        )


@requires(GenotypeHaplotypeCallerGvcf, CallVariantsWithHaplotypeCaller,
          FetchReferenceFasta, FetchEvaluationIntervalList,
          CreateSequenceDictionary)
class CNNScoreVariants(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        output_vcf_path = re.sub(
            r'\.vcf\.gz$', '.cnn.vcf.gz', self.input()[0][0].path
        )
        return [luigi.LocalTarget(output_vcf_path + s) for s in ['', '.tbi']]

    def run(self):
        output_cnn_vcf = Path(self.output()[0].path)
        run_id = '.'.join(output_cnn_vcf.name.split('.')[:-4])
        self.print_log(f'Score variants with CNN:\t{run_id}')
        gatk = self.cf['gatk']
        python3 = self.cf['python3']
        raw_vcf = Path(self.input()[0][0].path)
        cram = Path(self.input()[1][2].path)
        fa = Path(self.input()[2][0].path)
        evaluation_interval = Path(self.input()[3].path)
        self.setup_shell(
            run_id=run_id, commands=[gatk, python3], cwd=output_cnn_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} CNNScoreVariants'
                + f' --reference {fa}'
                + f' --input {cram}'
                + f' --variant {raw_vcf}'
                + f' --intervals {evaluation_interval}'
                + f' --output {output_cnn_vcf}'
                + ' --tensor-type read_tensor'
                + ' --disable-bam-index-caching '
                + str(self.cf['save_memory']).lower()
            ),
            input_files_or_dirs=[raw_vcf, fa, cram, evaluation_interval],
            output_files_or_dirs=[output_cnn_vcf, f'{output_cnn_vcf}.tbi']
        )


@requires(CNNScoreVariants, FetchHapmapVcf, FetchMillsIndelVcf,
          FetchEvaluationIntervalList)
class FilterVariantTranches(VclineTask):
    cf = luigi.DictParameter()
    snp_tranche = luigi.ListParameter(default=[99.9, 99.95])
    indel_tranche = luigi.ListParameter(default=[99.0, 99.4])
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        output_vcf_path = re.sub(
            r'\.vcf\.gz$', '.filtered.vcf.gz', self.input()[0][0].path
        )
        return [luigi.LocalTarget(output_vcf_path + s) for s in ['', '.tbi']]

    def run(self):
        output_filtered_vcf = Path(self.output()[0].path)
        run_id = '.'.join(output_filtered_vcf.name.split('.')[:-5])
        self.print_log(f'Apply tranche filtering:\t{run_id}')
        cnn_vcf = Path(self.input()[0][0].path)
        resource_vcfs = [
            Path(o.path) for o in [self.input()[1][0], self.input()[2][0]]
        ]
        evaluation_interval = Path(self.input()[3].path)
        gatk = self.cf['gatk']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=output_filtered_vcf.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} FilterVariantTranches'
                + f' --variant {cnn_vcf}'
                + ''.join([f' --resource {p}' for p in resource_vcfs])
                + f' --intervals {evaluation_interval}'
                + f' --output {output_filtered_vcf}'
                + ' --info-key CNN_2D'
                + ''.join(
                    [f' --snp-tranche {v}' for v in self.snp_tranche]
                    + [f' --indel-tranche {v}' for v in self.indel_tranche]
                )
                + ' --invalidate-previous-filters'
                + ' --disable-bam-index-caching '
                + str(self.cf['save_memory']).lower()
            ),
            input_files_or_dirs=[cnn_vcf, *resource_vcfs, evaluation_interval],
            output_files_or_dirs=[
                output_filtered_vcf, f'{output_filtered_vcf}.tbi'
            ]
        )


if __name__ == '__main__':
    luigi.run()
