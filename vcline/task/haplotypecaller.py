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
from .samtools import MergeSAMsIntoSortedSAM, SamtoolsViewAndSamtoolsIndex


@requires(FetchEvaluationIntervalList, FetchReferenceFASTA)
class SplitEvaluationIntervals(ShellTask):
    scatter_count = luigi.IntParameter(default=2)
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['ref_dir_path']).joinpath(
                    f'gatk/{i:04d}-scattered.interval_list'
                )
            ) for i in range(self.scatter_count)
        ]

    def run(self):
        input_interval_path = self.input()[0].path
        run_id = Path(input_interval_path).stem
        self.print_log(f'Split an evaluation interval list:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        fa_path = self.input()[1][0].path
        output_interval_paths = [o.path for o in self.output()]
        output_dir = Path(output_interval_paths[0]).parent
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=output_dir, remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} SplitIntervals'
                + f' --reference {fa_path}'
                + f' --intervals {input_interval_path}'
                + f' --scatter-count {self.scatter_count}'
                + f' --output {output_dir}'
            ),
            input_files_or_dirs=[input_interval_path, fa_path],
            output_files_or_dirs=output_interval_paths
        )


class PrepareEvaluationIntervals(luigi.WrapperTask):
    evaluation_interval_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def requires(self):
        return (
            SplitEvaluationIntervals(
                evaluation_interval_path=self.evaluation_interval_path,
                ref_fa_path=self.ref_fa_path,
                scatter_count=self.cf['n_worker'], cf=self.cf
            ) if self.cf['n_worker'] > 1 else [
                FetchEvaluationIntervalList(
                    evaluation_interval_path=self.evaluation_interval_path,
                    cf=self.cf
                )
            ]
        )

    def output(self):
        return self.input()


@requires(PrepareCRAMNormal, FetchReferenceFASTA, FetchDbsnpVCF,
          PrepareEvaluationIntervals, CreateSequenceDictionary)
class CallVariantsWithHaplotypeCaller(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['germline_snv_indel_gatk_dir_path']).joinpath(
                    Path(self.input()[0][0].path).stem
                    + f'.haplotypecaller.{s}'
                )
            ) for s in ['g.vcf.gz', 'g.vcf.gz.tbi', 'cram', 'cram.crai']
        ]

    def run(self):
        output_gvcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_gvcf_path).name.split('.')[:-4])
        self.print_log(
            f'Call germline variants with HaplotypeCaller:\t{run_id}'
        )
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        input_cram_path = self.input()[0][0].path
        fa_path = self.input()[1][0].path
        dbsnp_vcf_path = self.input()[2][0].path
        evaluation_interval_paths = [i.path for i in self.input()[3]]
        output_cram_path = self.output()[2].path
        output_path_prefix = '.'.join(output_gvcf_path.split('.')[:-3])
        if len(evaluation_interval_paths) == 1:
            tmp_prefixes = [output_path_prefix]
        else:
            tmp_prefixes = [
                '{0}.{1}'.format(output_path_prefix, Path(i).stem)
                for i in evaluation_interval_paths
            ]
        targets = yield [
            HaplotypeCaller(
                input_cram_path=input_cram_path, fa_path=fa_path,
                dbsnp_vcf_path=dbsnp_vcf_path, evaluation_interval_path=i,
                output_path_prefix=s, cf=self.cf
            ) for i, s in zip(evaluation_interval_paths, tmp_prefixes)
        ]
        if len(evaluation_interval_paths) == 1:
            yield SamtoolsViewAndSamtoolsIndex(
                input_sam_path=f'{tmp_prefixes[0]}.bam',
                output_sam_path=output_cram_path, fa_path=fa_path,
                samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'], remove_input=True,
                log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed']
            )
        else:
            self.setup_shell(
                run_id=run_id, log_dir_path=self.cf['log_dir_path'],
                commands=gatk, cwd=self.cf['germline_snv_indel_gatk_dir_path'],
                remove_if_failed=self.cf['remove_if_failed']
            )
            tmp_gvcf_paths = [f'{s}.g.vcf.gz' for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {gatk}{gatk_opts} CombineGVCFs'
                    + f' --reference {fa_path}'
                    + ''.join([f' --variant {p}' for p in tmp_gvcf_paths])
                    + f' --output {output_gvcf_path}'
                ),
                input_files_or_dirs=[*tmp_gvcf_paths, fa_path],
                output_files_or_dirs=[
                    output_gvcf_path, f'{output_gvcf_path}.tbi'
                ]
            )
            yield MergeSAMsIntoSortedSAM(
                input_sam_paths=[f'{s}.bam' for s in tmp_prefixes],
                output_sam_path=output_cram_path, fa_path=fa_path,
                samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'],
                memory_per_thread=self.cf['samtools_memory_per_thread'],
                index_sam=True, remove_input=False,
                log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed']
            )
            for t in targets:
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
        self.setup_shell(
            run_id='.'.join(Path(output_file_paths[0]).name.split('.')[:-3]),
            log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=Path(output_file_paths[0]).parent,
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} HaplotypeCaller'
                + f' --input {self.input_cram_path}'
                + f' --reference {self.fa_path}'
                + f' --dbsnp {self.dbsnp_vcf_path}'
                + f' --intervals {self.evaluation_interval_path}'
                + f' --output {output_file_paths[0]}'
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
        return [
            luigi.LocalTarget(
                re.sub(r'\.g\.vcf\.gz$', s, self.input()[0][0].path)
            ) for s in ['.vcf.gz', '.vcf.gz.tbi']
        ]

    def run(self):
        vcf_path = self.output()[0].path
        run_id = '.'.join(Path(vcf_path).name.split('.')[:-3])
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
            cwd=self.cf['germline_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} GenotypeGVCFs'
                + f' --reference {fa_path}'
                + f' --variant {gvcf_path}'
                + f' --dbsnp {dbsnp_vcf_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {vcf_path}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                gvcf_path, fa_path, dbsnp_vcf_path, evaluation_interval_path
            ],
            output_files_or_dirs=[vcf_path, f'{vcf_path}.tbi']
        )


@requires(GenotypeHaplotypeCallerGVCF, CallVariantsWithHaplotypeCaller,
          FetchReferenceFASTA, FetchEvaluationIntervalList,
          CreateSequenceDictionary)
class CNNScoreVariants(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                re.sub(r'\.vcf\.gz$', f'.cnn.{s}', self.input()[0][0].path)
            ) for s in ['vcf.gz', 'vcf.gz.tbi']
        ]

    def run(self):
        cnn_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(cnn_vcf_path).name.split('.')[:-4])
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
            commands=[gatk, python3],
            cwd=self.cf['germline_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CNNScoreVariants'
                + f' --reference {fa_path}'
                + f' --input {cram_path}'
                + f' --variant {raw_vcf_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {cnn_vcf_path}'
                + ' --tensor-type read_tensor'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                raw_vcf_path, fa_path, cram_path, evaluation_interval_path
            ],
            output_files_or_dirs=[cnn_vcf_path, f'{cnn_vcf_path}.tbi']
        )


@requires(CNNScoreVariants, FetchHapmapVCF, FetchMillsIndelVCF,
          FetchEvaluationIntervalList)
class FilterVariantTranches(ShellTask):
    cf = luigi.DictParameter()
    snp_tranche = luigi.ListParameter(default=[99.9, 99.95])
    indel_tranche = luigi.ListParameter(default=[99.0, 99.4])
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                re.sub(
                    r'\.vcf\.gz$', f'.filtered.{s}', self.input()[0][0].path
                )
            ) for s in ['vcf.gz', 'vcf.gz.tbi']
        ]

    def run(self):
        filtered_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(filtered_vcf_path).name.split('.')[:-5])
        self.print_log(f'Apply tranche filtering:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        cnn_vcf_path = self.input()[0][0].path
        resource_vcf_paths = [self.input()[1][0].path, self.input()[2][0].path]
        evaluation_interval_path = self.input()[3].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['germline_snv_indel_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} FilterVariantTranches'
                + f' --variant {cnn_vcf_path}'
                + ''.join([f' --resource {p}' for p in resource_vcf_paths])
                + f' --intervals {evaluation_interval_path}'
                + f' --output {filtered_vcf_path}'
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
                filtered_vcf_path, f'{filtered_vcf_path}.tbi'
            ]
        )


if __name__ == '__main__':
    luigi.run()
