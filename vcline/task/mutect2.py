#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .base import ShellTask
from .ref import (CreateGnomadBiallelicSnpVCF, CreateSequenceDictionary,
                  FetchEvaluationIntervalList, FetchGnomadVCF,
                  FetchReferenceFASTA, PrepareEvaluationIntervals)
from .samtools import MergeSAMsIntoSortedSAM, SamtoolsViewAndSamtoolsIndex


class GetPileupSummaries(ShellTask):
    cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    gnomad_common_biallelic_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['somatic_snv_indel_gatk_dir_path']).joinpath(
                Path(self.cram_path).stem + '.pileup.table'
            )
        )

    def run(self):
        run_id = Path(self.cram_path).stem
        self.print_log(f'Get pileup summary:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        pileup_table_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} GetPileupSummaries'
                + f' --input {self.cram_path}'
                + f' --reference {self.fa_path}'
                + f' --variant {self.gnomad_common_biallelic_vcf_path}'
                + f' --intervals {self.evaluation_interval_path}'
                + f' --output {pileup_table_path}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                self.cram_path, self.fa_path, self.evaluation_interval_path,
                self.gnomad_common_biallelic_vcf_path
            ],
            output_files_or_dirs=pileup_table_path
        )


@requires(PrepareCRAMTumor, PrepareCRAMNormal, FetchReferenceFASTA,
          FetchEvaluationIntervalList, CreateGnomadBiallelicSnpVCF,
          CreateSequenceDictionary)
class CalculateContamination(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['somatic_snv_indel_gatk_dir_path']).joinpath(
                    create_matched_id(*[i[0].path for i in self.input()[0:2]])
                    + s
                )
            ) for s in ['.contamination.table', '.segment.table']
        ]

    def run(self):
        input_targets = yield [
            GetPileupSummaries(
                cram_path=self.input()[i][0].path,
                fa_path=self.input()[2][0].path,
                evaluation_interval_path=self.input()[3].path,
                gnomad_common_biallelic_vcf_path=self.input()[4][0].path,
                cf=self.cf
            ) for i in range(2)
        ]
        contamination_table_path = self.output()[0].path
        run_id = '.'.join(Path(contamination_table_path).name.split('.')[:-2])
        self.print_log(f'Calculate cross-sample contamination:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        pileup_table_paths = [i.path for i in input_targets]
        segment_table_path = self.output()[1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
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


@requires(PrepareCRAMTumor, PrepareCRAMNormal, FetchReferenceFASTA,
          PrepareEvaluationIntervals, FetchGnomadVCF,
          CreateSequenceDictionary)
class CallVariantsWithMutect2(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['somatic_snv_indel_gatk_dir_path']).joinpath(
                    create_matched_id(*[i[0].path for i in self.input()[0:2]])
                    + f'.mutect2.{s}'
                )
            ) for s in [
                'vcf.gz', 'vcf.gz.tbi', 'vcf.gz.stats', 'cram', 'cram.crai',
                'read-orientation-model.tar.gz'
            ]
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Call somatic variants with Mutect2:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        input_cram_paths = [i[0].path for i in self.input()[0:2]]
        fa_path = self.input()[2][0].path
        evaluation_interval_paths = [i.path for i in self.input()[3]]
        gnomad_vcf_path = self.input()[4][0].path
        output_stats_path = self.output()[2].path
        output_cram_path = self.output()[3].path
        ob_priors_path = self.output()[5].path
        output_path_prefix = '.'.join(output_vcf_path.split('.')[:-2])
        if len(evaluation_interval_paths) == 1:
            tmp_prefixes = [output_path_prefix]
        else:
            tmp_prefixes = [
                '{0}.{1}'.format(output_path_prefix, Path(i).stem)
                for i in evaluation_interval_paths
            ]
        targets = yield [
            Mutect2(
                input_cram_paths=input_cram_paths, fa_path=fa_path,
                gnomad_vcf_path=gnomad_vcf_path, evaluation_interval_path=i,
                normal_name=self.sample_names[1], output_path_prefix=s,
                cf=self.cf
            ) for i, s in zip(evaluation_interval_paths, tmp_prefixes)
        ]
        f1r2_paths = [f'{s}.f1r2.tar.gz' for s in tmp_prefixes]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} LearnReadOrientationModel'
                + ''.join([f' --input {f}' for f in f1r2_paths])
                + f' --output {ob_priors_path}'
            ),
            input_files_or_dirs=f1r2_paths, output_files_or_dirs=ob_priors_path
        )
        if len(evaluation_interval_paths) == 1:
            yield SamtoolsViewAndSamtoolsIndex(
                input_sam_path=f'{tmp_prefixes[0]}.bam',
                output_sam_path=output_cram_path, fa_path=fa_path,
                samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'], remove_input=True,
                log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed'],
                quiet=self.cf['quiet']
            )
        else:
            tmp_vcf_paths = [f'{s}.vcf.gz' for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {gatk}{gatk_opts} MergeVcfs'
                    + ''.join([f' --INPUT {v}' for v in tmp_vcf_paths])
                    + f' --OUTPUT {output_vcf_path}'
                ),
                input_files_or_dirs=tmp_vcf_paths,
                output_files_or_dirs=output_vcf_path
            )
            tmp_stats_paths = [f'{s}.vcf.gz.stats' for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {gatk}{gatk_opts} MergeMutectStats'
                    + ''.join([f' --stats {s}' for s in tmp_stats_paths])
                    + f' --output {output_stats_path}'
                ),
                input_files_or_dirs=tmp_stats_paths,
                output_files_or_dirs=output_stats_path
            )
            yield MergeSAMsIntoSortedSAM(
                input_sam_paths=[f'{s}.bam' for s in tmp_prefixes],
                output_sam_path=output_cram_path, fa_path=fa_path,
                samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'],
                memory_per_thread=self.cf['samtools_memory_per_thread'],
                index_sam=True, remove_input=False,
                log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed'],
                quiet=self.cf['quiet']
            )
            for t in targets:
                tmp_file_paths = [o.path for o in t]
                self.run_shell(
                    args=('rm -f ' + ' '.join(tmp_file_paths)),
                    input_files_or_dirs=tmp_file_paths
                )


class Mutect2(ShellTask):
    input_cram_paths = luigi.ListParameter()
    fa_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    normal_name = luigi.Parameter()
    output_path_prefix = luigi.Parameter()
    cf = luigi.DictParameter()
    message = luigi.Parameter(default='')
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
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        n_cpu = self.cf['n_cpu_per_worker']
        output_file_paths = [o.path for o in self.output()]
        self.setup_shell(
            run_id='.'.join(Path(output_file_paths[0]).name.split('.')[:-2]),
            log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=Path(output_file_paths[0]).parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} Mutect2'
                + ''.join([f' --input {p}' for p in self.input_cram_paths])
                + f' --reference {self.fa_path}'
                + f' --intervals {self.evaluation_interval_path}'
                + f' --germline-resource {self.gnomad_vcf_path}'
                + f' --output {output_file_paths[0]}'
                + f' --bam-output {output_file_paths[3]}'
                + f' --f1r2-tar-gz {output_file_paths[4]}'
                + f' --normal-sample {self.normal_name}'
                + ' --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP'
                + f' --native-pair-hmm-threads {n_cpu}'
                + ' --max-mnp-distance 0'
                + ' --create-output-bam-index false'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                *self.input_cram_paths, self.fa_path,
                self.evaluation_interval_path, self.gnomad_vcf_path
            ],
            output_files_or_dirs=output_file_paths
        )


@requires(CallVariantsWithMutect2, FetchReferenceFASTA,
          FetchEvaluationIntervalList, CalculateContamination,
          CreateSequenceDictionary)
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
        fa_path = self.input()[1][0].path
        evaluation_interval_path = self.input()[2].path
        filtering_stats_path = self.output()[2].path
        contamination_table_path = self.input()[3][0].path
        segment_table_path = self.input()[3][1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
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
