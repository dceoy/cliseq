#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareNormalCRAM, PrepareTumorCRAM
from .base import ShellTask
from .haplotypecaller import PrepareEvaluationIntervals
from .ref import (CreateFASTAIndex, CreateGnomadBiallelicSnpVCF,
                  FetchEvaluationIntervalList, FetchGnomadVCF,
                  FetchReferenceFASTA)
from .samtools import BAM2CRAM, MergeBAMsIntoCRAM, SamtoolsIndex


@requires(FetchReferenceFASTA, FetchEvaluationIntervalList,
          CreateGnomadBiallelicSnpVCF)
class GetPileupSummaries(ShellTask):
    cram_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['somatic_snv_indel_gatk_dir_path']).joinpath(
                    Path(self.cram_path).stem + '.pileup.table'
                )
            )
        )

    def run(self):
        run_id = Path(self.cram_path).stem
        self.print_log(f'Get pileup summary:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        pileup_table_path = self.output().path
        fa_path = self.input()[0].path
        evaluation_interval_path = self.input()[1].path
        gnomad_common_biallelic_vcf_path = self.input()[2][0].path
        pileup_table_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} GetPileupSummaries'
                + f' --reference {fa_path}'
                + f' --input {self.cram_path}'
                + f' --variant {gnomad_common_biallelic_vcf_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {pileup_table_path}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                self.cram_path, fa_path, evaluation_interval_path,
                gnomad_common_biallelic_vcf_path
            ],
            output_files_or_dirs=pileup_table_path
        )


@requires(PrepareTumorCRAM, PrepareNormalCRAM)
class GetPileupSummariesMatched(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            GetPileupSummaries(
                cram_path=i[0].path, ref_fa_paths=self.ref_fa_paths,
                gnomad_vcf_path=self.gnomad_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            ).output() for i in self.input()
        ]


@requires(GetPileupSummariesMatched)
class CalculateContamination(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['somatic_snv_indel_gatk_dir_path']).joinpath(
                        create_matched_id(*[i.path for i in self.input()]) + s
                    )
                )
            ) for s in ['.contamination.table', '.segment.table']
        ]

    def run(self):
        contamination_table_path = self.output()[0].path
        run_id = '.'.join(Path(contamination_table_path).name.split('.')[:-2])
        self.print_log(f'Calculate cross-sample contamination:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        pileup_table_paths = [i.path for i in self.input()]
        segment_table_path = self.output()[1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CalculateContamination'
                + f' --input {pileup_table_paths[0]}'
                + f' --matched-normal {pileup_table_paths[1]}'
                + f' --output {contamination_table_path}'
                + f' --tumor-segmentation {segment_table_path}'
            ),
            input_files_or_dirs=pileup_table_paths,
            output_files_or_dirs=[contamination_table_path, segment_table_path]
        )


@requires(PrepareTumorCRAM, PrepareNormalCRAM, FetchReferenceFASTA,
          CreateFASTAIndex, PrepareEvaluationIntervals, FetchGnomadVCF)
class CallVariantsWithMutect2(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['somatic_snv_indel_gatk_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0:2]]
                        ) + f'.mutect2.{s}'
                    )
                )
            ) for s in [
                'vcf.gz', 'vcf.gz.tbi', 'vcf.gz.stats', 'cram', 'cram.crai',
                'read-orientation-model.tar.gz'
            ]
        ]

    def run(self):
        raw_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(raw_vcf_path).name.split('.')[:-3])
        self.print_log(f'Call somatic variants with Mutect2:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        samtools = self.cf['samtools']
        save_memory = str(self.cf['save_memory']).lower()
        n_cpu = self.cf['n_cpu_per_worker']
        memory_per_thread = self.cf['samtools_memory_per_thread']
        input_cram_paths = [i[0].path for i in self.input()[0:2]]
        fa_path = self.input()[2].path
        evaluation_interval_paths = [i.path for i in self.input()[4]]
        gnomad_vcf_path = self.input()[5][0].path
        raw_stats_path = self.output()[2].path
        output_cram_path = self.output()[3].path
        ob_priors_path = self.output()[5].path
        if len(evaluation_interval_paths) == 1:
            tmp_bam_paths = [re.sub(r'(\.cram)$', '.bam', output_cram_path)]
            tmp_vcf_paths = [raw_vcf_path]
        else:
            tmp_bam_paths = [
                re.sub(
                    r'\.cram$', '.{}.bam'.format(Path(i).stem),
                    output_cram_path
                ) for i in evaluation_interval_paths
            ]
            tmp_vcf_paths = [
                re.sub(
                    r'\.vcf\.gz$', '.{}.vcf.gz'.format(Path(i).stem),
                    raw_vcf_path
                ) for i in evaluation_interval_paths
            ]
        tmp_tbi_paths = [f'{p}.tbi' for p in tmp_vcf_paths]
        f1r2_paths = [
            re.sub(
                r'\.cram$', '.{}.f1r2.tar.gz'.format(Path(i).stem),
                output_cram_path
            ) for i in evaluation_interval_paths
        ]
        tmp_stats_paths = [f'{v}.stats' for v in tmp_vcf_paths]
        normal_name = self.sample_names[1]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, samtools],
            cwd=self.cf['somatic_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {gatk}{gatk_opts} Mutect2'
                    + f' --reference {fa_path}'
                    + ''.join([f' --input {p}' for p in input_cram_paths])
                    + f' --intervals {i}'
                    + f' --germline-resource {gnomad_vcf_path}'
                    + f' --output {v}'
                    + f' --bam-output {b}'
                    + f' --f1r2-tar-gz {f}'
                    + f' --normal-sample {normal_name}'
                    + ' --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP'
                    + f' --native-pair-hmm-threads {n_cpu}'
                    + f' --disable-bam-index-caching {save_memory}'
                    + ' --max-mnp-distance 0'
                    + ' --create-output-bam-index false'
                ) for i, v, b, f in zip(
                    evaluation_interval_paths, tmp_vcf_paths, tmp_bam_paths,
                    f1r2_paths
                )
            ],
            input_files_or_dirs=[
                *input_cram_paths, fa_path, *evaluation_interval_paths,
                gnomad_vcf_path
            ],
            output_files_or_dirs=[
                *tmp_vcf_paths, *tmp_tbi_paths, *tmp_bam_paths, *f1r2_paths,
                *tmp_stats_paths
            ],
            asynchronous=(len(evaluation_interval_paths) > 1)
        )
        if len(tmp_bam_paths) > 1:
            yield MergeBAMsIntoCRAM(
                bam_paths=tmp_bam_paths, output_cram_path=output_cram_path,
                fa_path=fa_path, samtools=samtools, n_cpu=n_cpu,
                memory_per_thread=memory_per_thread,
                log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed']
            )
        else:
            yield BAM2CRAM(
                sam_path=tmp_bam_paths[0], fa_path=fa_path, samtools=samtools,
                n_cpu=n_cpu, log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed']
            )
        yield SamtoolsIndex(
            sam_path=output_cram_path, samtools=samtools, n_cpu=n_cpu,
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} LearnReadOrientationModel'
                + ''.join([f' --input {f}' for f in f1r2_paths])
                + f' --output {ob_priors_path}'
            ),
            input_files_or_dirs=f1r2_paths, output_files_or_dirs=ob_priors_path
        )
        if len(tmp_vcf_paths) > 1:
            self.run_shell(
                args=[
                    (
                        f'set -e && {gatk}{gatk_opts} MergeMutectStats'
                        + ''.join([f' --stats {s}' for s in tmp_stats_paths])
                        + f' --output {raw_stats_path}'
                    ),
                    ('rm -f ' + ' '.join(tmp_stats_paths))
                ],
                input_files_or_dirs=tmp_stats_paths,
                output_files_or_dirs=raw_stats_path
            )
            self.run_shell(
                args=[
                    (
                        f'set -e && {gatk}{gatk_opts} MergeVcfs'
                        + ''.join([f' --INPUT {v}' for v in tmp_vcf_paths])
                        + f' --OUTPUT {raw_vcf_path}'
                    ),
                    ('rm -f ' + ' '.join([*tmp_vcf_paths, *tmp_tbi_paths]))
                ],
                input_files_or_dirs=[*tmp_vcf_paths, *tmp_tbi_paths],
                output_files_or_dirs=raw_vcf_path
            )


@requires(CallVariantsWithMutect2, FetchReferenceFASTA,
          FetchEvaluationIntervalList, CalculateContamination)
class FilterMutectCalls(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                re.sub(
                    r'\.vcf\.gz$', f'.filtered.{s}', self.input()[0][0].path
                )
            ) for s in ['vcf.gz', 'vcf.gz.tbi', 'vcf.gz.stats']
        ]

    def run(self):
        filtered_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(filtered_vcf_path).name.split('.')[:-4])
        self.print_log(f'Filter somatic variants called by Mutect2:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        raw_vcf_path = self.input()[0][0].path
        raw_stats_path = self.input()[0][2].path
        ob_priors_path = self.input()[0][5].path
        fa_path = self.input()[1].path
        evaluation_interval_path = self.input()[2].path
        filtering_stats_path = self.output()[2].path
        contamination_table_path = self.input()[3][0].path
        segment_table_path = self.input()[3][1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} FilterMutectCalls'
                + f' --reference {fa_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --variant {raw_vcf_path}'
                + f' --stats {raw_stats_path}'
                + f' --contamination-table {contamination_table_path}'
                + f' --tumor-segmentation {segment_table_path}'
                + f' --orientation-bias-artifact-priors {ob_priors_path}'
                + f' --output {filtered_vcf_path}'
                + f' --filtering-stats {filtering_stats_path}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                raw_vcf_path, fa_path, evaluation_interval_path,
                raw_stats_path, ob_priors_path,
                contamination_table_path, segment_table_path,
            ],
            output_files_or_dirs=[filtered_vcf_path, filtering_stats_path]
        )


if __name__ == '__main__':
    luigi.run()
