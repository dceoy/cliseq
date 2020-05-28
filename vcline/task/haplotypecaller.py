#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from .align import PrepareCRAMNormal
from .base import ShellTask
from .ref import (CreateSequenceDictionary, FetchDbsnpVCF,
                  FetchEvaluationIntervalList, FetchHapmapVCF,
                  FetchMillsIndelVCF, FetchReferenceFASTA)
from .samtools import samtools_merge_and_index, samtools_view_and_index


class SplitIntervals(ShellTask):
    interval_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter()
    scatter_count = luigi.IntParameter(default=2)
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        if self.scatter_count > 1:
            return [
                luigi.LocalTarget(
                    Path(self.dest_dir_path).joinpath(
                        f'{i:04d}-scattered.interval_list'
                    )
                ) for i in range(self.scatter_count)
            ]
        else:
            return [luigi.LocalTarget(self.interval_path)]

    def run(self):
        run_id = Path(self.interval_path).stem
        self.print_log(f'Split an interval list:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        output_interval_paths = [o.path for o in self.output()]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.dest_dir_path,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} SplitIntervals'
                + f' --reference {self.fa_path}'
                + f' --intervals {self.interval_path}'
                + f' --scatter-count {self.scatter_count}'
                + f' --output {self.dest_dir_path}'
            ),
            input_files_or_dirs=[self.interval_path, self.fa_path],
            output_files_or_dirs=output_interval_paths
        )


@requires(PrepareCRAMNormal, FetchReferenceFASTA, FetchDbsnpVCF,
          FetchEvaluationIntervalList, CreateSequenceDictionary)
class CallVariantsWithHaplotypeCaller(ShellTask):
    cf = luigi.DictParameter()
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
        scatter_count = self.cf['n_worker']
        fa_path = self.input()[1][0].path
        interval_targets = yield SplitIntervals(
            interval_path=self.input()[3].path, fa_path=fa_path,
            dest_dir_path=str(output_gvcf.parent.parent),
            scatter_count=scatter_count, cf=self.cf
        )
        input_cram_path = self.input()[0][0].path
        dbsnp_vcf_path = self.input()[2][0].path
        output_path_prefix = '.'.join(str(output_gvcf).split('.')[:-3])
        if scatter_count == 1:
            tmp_prefixes = [output_path_prefix]
        else:
            tmp_prefixes = [
                '{0}.{1}'.format(output_path_prefix, Path(i.path).stem)
                for i in interval_targets
            ]
        input_targets = yield [
            HaplotypeCaller(
                input_cram_path=input_cram_path, fa_path=fa_path,
                dbsnp_vcf_path=dbsnp_vcf_path, evaluation_interval_path=i.path,
                output_path_prefix=s, cf=self.cf
            ) for i, s in zip(interval_targets, tmp_prefixes)
        ]
        run_id = '.'.join(output_gvcf.name.split('.')[:-4])
        self.print_log(
            f'Call germline variants with HaplotypeCaller:\t{run_id}'
        )
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        samtools = self.cf['samtools']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_mb = self.cf['memory_mb_per_worker']
        output_cram_path = self.output()[2].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=output_gvcf.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        if scatter_count == 1:
            tmp_bam_path = f'{tmp_prefixes[0]}.bam'
            samtools_view_and_index(
                shelltask=self, samtools=samtools, input_sam_path=tmp_bam_path,
                fa_path=fa_path, output_sam_path=output_cram_path, n_cpu=n_cpu
            )
            self.run_shell(
                args=f'rm -f {tmp_bam_path}', input_files_or_dirs=tmp_bam_path
            )
        else:
            tmp_gvcf_paths = [f'{s}.g.vcf.gz' for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {gatk}{gatk_opts} CombineGVCFs'
                    + f' --reference {fa_path}'
                    + ''.join([f' --variant {p}' for p in tmp_gvcf_paths])
                    + f' --output {output_gvcf}'
                ),
                input_files_or_dirs=[*tmp_gvcf_paths, fa_path],
                output_files_or_dirs=[output_gvcf, f'{output_gvcf}.tbi']
            )
            samtools_merge_and_index(
                shelltask=self, samtools=samtools,
                input_sam_paths=[f'{s}.bam' for s in tmp_prefixes],
                fa_path=fa_path, output_sam_path=output_cram_path, n_cpu=n_cpu,
                memory_mb=memory_mb
            )
            for t in input_targets:
                tmp_file_paths = [o.path for o in t]
                self.run_shell(
                    args=('rm -f ' + ' '.join(tmp_file_paths)),
                    input_files_or_dirs=tmp_file_paths
                )


class HaplotypeCaller(ShellTask):
    input_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dbsnp_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    output_path_prefix = luigi.Parameter()
    cf = luigi.DictParameter()
    message = luigi.Parameter(default='')
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(f'{self.output_path_prefix}.{s}')
            for s in ['g.vcf.gz', 'g.vcf.gz.tbi', 'bam']
        ]

    def run(self):
        if self.message:
            self.print_log(self.message)
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        n_cpu = self.cf['n_cpu_per_worker']
        output_file_paths = [o.path for o in self.output()]
        output_gvcf = Path(output_file_paths[0])
        self.setup_shell(
            run_id='.'.join(output_gvcf.name.split('.')[:-3]),
            log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=output_gvcf.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} HaplotypeCaller'
                + f' --input {self.input_cram_path}'
                + f' --reference {self.fa_path}'
                + f' --dbsnp {self.dbsnp_vcf_path}'
                + f' --intervals {self.evaluation_interval_path}'
                + f' --output {output_gvcf}'
                + f' --bam-output {output_file_paths[2]}'
                + ' --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP'
                + f' --native-pair-hmm-threads {n_cpu}'
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
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                self.input_cram_path, self.fa_path, self.dbsnp_vcf_path,
                self.evaluation_interval_path
            ],
            output_files_or_dirs=output_file_paths
        )


@requires(CallVariantsWithHaplotypeCaller, FetchReferenceFASTA,
          FetchDbsnpVCF, FetchEvaluationIntervalList, CreateSequenceDictionary)
class GenotypeHaplotypeCallerGVCF(ShellTask):
    cf = luigi.DictParameter()
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
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        gvcf_path = self.input()[0][0].path
        fa_path = self.input()[1][0].path
        dbsnp_vcf_path = self.input()[2][0].path
        evaluation_interval_path = self.input()[3].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=output_vcf.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} GenotypeGVCFs'
                + f' --reference {fa_path}'
                + f' --variant {gvcf_path}'
                + f' --dbsnp {dbsnp_vcf_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {output_vcf}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                gvcf_path, fa_path, dbsnp_vcf_path, evaluation_interval_path
            ],
            output_files_or_dirs=[output_vcf, f'{output_vcf}.tbi']
        )


@requires(GenotypeHaplotypeCallerGVCF, CallVariantsWithHaplotypeCaller,
          FetchReferenceFASTA, FetchEvaluationIntervalList,
          CreateSequenceDictionary)
class CNNScoreVariants(ShellTask):
    cf = luigi.DictParameter()
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
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        python3 = self.cf['python3']
        save_memory = str(self.cf['save_memory']).lower()
        raw_vcf_path = self.input()[0][0].path
        cram_path = self.input()[1][2].path
        fa_path = self.input()[2][0].path
        evaluation_interval_path = self.input()[3].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, python3], cwd=output_cnn_vcf.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CNNScoreVariants'
                + f' --reference {fa_path}'
                + f' --input {cram_path}'
                + f' --variant {raw_vcf_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {output_cnn_vcf}'
                + ' --tensor-type read_tensor'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                raw_vcf_path, fa_path, cram_path, evaluation_interval_path
            ],
            output_files_or_dirs=[output_cnn_vcf, f'{output_cnn_vcf}.tbi']
        )


@requires(CNNScoreVariants, FetchHapmapVCF, FetchMillsIndelVCF,
          FetchEvaluationIntervalList)
class FilterVariantTranches(ShellTask):
    cf = luigi.DictParameter()
    snp_tranche = luigi.ListParameter(default=[99.9, 99.95])
    indel_tranche = luigi.ListParameter(default=[99.0, 99.4])
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
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        cnn_vcf_path = self.input()[0][0].path
        resource_vcf_paths = [self.input()[1][0].path, self.input()[2][0].path]
        evaluation_interval_path = self.input()[3].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=output_filtered_vcf.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} FilterVariantTranches'
                + f' --variant {cnn_vcf_path}'
                + ''.join([f' --resource {p}' for p in resource_vcf_paths])
                + f' --intervals {evaluation_interval_path}'
                + f' --output {output_filtered_vcf}'
                + ' --info-key CNN_2D'
                + ''.join(
                    [f' --snp-tranche {v}' for v in self.snp_tranche]
                    + [f' --indel-tranche {v}' for v in self.indel_tranche]
                )
                + ' --invalidate-previous-filters'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                cnn_vcf_path, *resource_vcf_paths, evaluation_interval_path
            ],
            output_files_or_dirs=[
                output_filtered_vcf, f'{output_filtered_vcf}.tbi'
            ]
        )


if __name__ == '__main__':
    luigi.run()
