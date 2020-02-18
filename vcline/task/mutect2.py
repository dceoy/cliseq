#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id, print_log
from .align import PrepareCRAMs
from .base import ShellTask
from .haplotypecaller import ApplyVQSR, PrepareEvaluationIntervals
from .ref import (CreateFASTAIndex, CreateGnomadBiallelicSnpVCF,
                  FetchGnomadVCF, FetchReferenceFASTA)


@requires(PrepareCRAMs, FetchReferenceFASTA, PrepareEvaluationIntervals,
          CreateGnomadBiallelicSnpVCF)
class CalculateContamination(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['mutect2_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0]]
                        ) + f'.{s}.table'
                    )
                )
            ) for s in ['contamination', 'segment']
        ]

    def run(self):
        contamination_table_path = self.output()[0].path
        run_id = '.'.join(Path(contamination_table_path).name.split('.')[:-2])
        print_log(f'Calculate cross-sample contamination:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        n_cpu = self.cf['n_cpu_per_worker']
        input_cram_paths = [i[0].path for i in self.input()[0]]
        fa_path = self.input()[1].path
        evaluation_interval_path = self.input()[2].path
        gnomad_common_biallelic_vcf_path = self.input()[3][0].path
        pileup_table_paths = [f'{p}.pileup.table' for p in input_cram_paths]
        segment_table_path = self.output()[1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['mutect2_dir_path']
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {gatk}{gatk_opts} GetPileupSummaries'
                    + f' --reference {fa_path}'
                    + f' --input {c}'
                    + f' --variant {gnomad_common_biallelic_vcf_path}'
                    + f' --intervals {evaluation_interval_path}'
                    + f' --output {t}'
                ) for c, t in zip(input_cram_paths, pileup_table_paths)
            ],
            input_files=[
                *input_cram_paths, fa_path, evaluation_interval_path,
                gnomad_common_biallelic_vcf_path
            ],
            output_files=pileup_table_paths, asynchronous=(n_cpu > 1)
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CalculateContamination'
                + f' --input {pileup_table_paths[0]}'
                + f' --matched-normal {pileup_table_paths[1]}'
                + f' --output {contamination_table_path}'
                + f' --tumor-segmentation {segment_table_path}'
            ),
            input_files=pileup_table_paths,
            output_files=[contamination_table_path, segment_table_path]
        )


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          PrepareEvaluationIntervals, FetchGnomadVCF)
class CallVariantsWithMutect2(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['mutect2_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0]]
                        ) + f'.Mutect2.{s}'
                    )
                )
            ) for s in [
                'raw.vcf.gz', 'raw.vcf.stats', 'cram', 'cram.crai',
                'read-orientation-model.tar.gz'
            ]
        ]

    def run(self):
        raw_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(raw_vcf_path).name.split('.')[:-4])
        print_log(f'Call somatic variants with Mutect2:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        samtools = self.cf['samtools']
        save_memory = str(self.cf['memory_mb_per_worker'] < 8 * 1024).lower()
        n_cpu = self.cf['n_cpu_per_worker']
        memory_per_thread = self.cf['samtools_memory_per_thread']
        input_cram_paths = [i[0].path for i in self.input()[0]]
        fa_path = self.input()[1].path
        fai_path = self.input()[2].path
        evaluation_interval_paths = [i.path for i in self.input()[3]]
        gnomad_vcf_path = self.input()[4][0].path
        raw_stats_path = self.output()[1].path
        output_cram_path = self.output()[2].path
        ob_priors_path = self.output()[4].path
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
                    r'\.raw\.vcf\.gz$', '.{}.raw.vcf.gz'.format(Path(i).stem),
                    raw_vcf_path
                ) for i in evaluation_interval_paths
            ]
        f1r2_paths = [
            re.sub(
                r'\.cram$', '.{}.f1r2.tar.gz'.format(Path(i).stem),
                output_cram_path
            ) for i in evaluation_interval_paths
        ]
        tmp_stats_paths = [
            re.sub(r'\.gz$', '.stats', v) for v in tmp_vcf_paths
        ]
        normal_name = self.sample_names[1]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, samtools], cwd=self.cf['mutect2_dir_path'],
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
            input_files=[
                *input_cram_paths, fa_path, *evaluation_interval_paths,
                gnomad_vcf_path
            ],
            output_files=[
                *tmp_vcf_paths, *tmp_bam_paths, *f1r2_paths, *tmp_stats_paths
            ],
            asynchronous=(len(evaluation_interval_paths) > 1)
        )
        self.run_shell(
            args=[
                (
                    (
                        'set -eo pipefail && '
                        + f'{samtools} merge -@ {n_cpu} -rh - '
                        + ' '.join(tmp_bam_paths)
                        + f' | {samtools} sort -@ {n_cpu}'
                        + f' -m {memory_per_thread}'
                        + f' -T {output_cram_path}.sort -'
                        + f' | {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                        + f' -o {output_cram_path} -'
                    ) if len(tmp_bam_paths) > 1 else (
                        'set -e && '
                        + f'{samtools} view -@ {n_cpu} -T {fa_path} -CS'
                        + f' -o {output_cram_path} {tmp_bam_paths[0]}'
                    )
                ),
                (f'set -e && rm -f ' + ' '.join(tmp_bam_paths))
            ],
            input_files=[*tmp_bam_paths, fa_path, fai_path],
            output_files=output_cram_path
        )
        self.run_shell(
            args=[
                f'set -e && {samtools} quickcheck -v {output_cram_path}',
                f'set -e && {samtools} index -@ {n_cpu} {output_cram_path}'
            ],
            input_files=output_cram_path,
            output_files=f'{output_cram_path}.crai'
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} LearnReadOrientationModel'
                + ''.join([f' --input {f}' for f in f1r2_paths])
                + f' --output {ob_priors_path}'
            ),
            input_files=f1r2_paths, output_files=ob_priors_path
        )
        if len(tmp_vcf_paths) > 1:
            self.run_shell(
                args=[
                    (
                        f'set -e && {gatk}{gatk_opts} MergeMutectStats'
                        + ''.join([f' --stats {s}' for s in tmp_stats_paths])
                        + f' --output {raw_stats_path}'
                    ),
                    (f'set -e && rm -f ' + ' '.join(tmp_stats_paths))
                ],
                input_files=tmp_stats_paths, output_files=raw_stats_path
            )
            self.run_shell(
                args=[
                    (
                        f'set -e && {gatk}{gatk_opts} MergeVcfs'
                        + ''.join([f' --INPUT {v}' for v in tmp_vcf_paths])
                        + f' --OUTPUT {raw_vcf_path}'
                    ),
                    (f'set -e && rm -f ' + ' '.join(tmp_vcf_paths))
                ],
                input_files=tmp_vcf_paths, output_files=raw_vcf_path
            )


@requires(CallVariantsWithMutect2, FetchReferenceFASTA,
          PrepareEvaluationIntervals, CalculateContamination)
class FilterMutect2Calls(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.raw\.vcf\.gz$', '.vcf', self.input()[0][0].path)
        )

    def run(self):
        filtered_vcf_path = self.output().path
        run_id = '.'.join(Path(filtered_vcf_path).name.split('.')[:-2])
        print_log(f'Filter somatic variants called by Mutect2:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        raw_vcf_path = self.input()[0][0].path
        raw_stats_path = self.input()[0][1].path
        ob_priors_path = self.input()[0][4].path
        fa_path = self.input()[1].path
        evaluation_interval_path = self.input()[2].path
        filtering_stats_path = re.sub(r'\.gz$', '.stats', filtered_vcf_path)
        contamination_table_path = self.input()[3][0].path
        segment_table_path = self.input()[3][1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['mutect2_dir_path']
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
            ),
            input_files=[
                raw_vcf_path, fa_path, evaluation_interval_path,
                raw_stats_path, ob_priors_path,
                contamination_table_path, segment_table_path,
            ],
            output_files=filtered_vcf_path
        )


class CallVariantsWithGATK(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    known_indel_vcf_paths = luigi.ListParameter()
    hapmap_vcf_path = luigi.Parameter()
    omni_vcf_path = luigi.Parameter()
    snp_1000g_vcf_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        return [
            ApplyVQSR(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                hapmap_vcf_path=self.hapmap_vcf_path,
                omni_vcf_path=self.omni_vcf_path,
                snp_1000g_vcf_path=self.snp_1000g_vcf_path, cf=self.cf
            ),
            FilterMutect2Calls(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                gnomad_vcf_path=self.gnomad_vcf_path, cf=self.cf
            )
        ]

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
