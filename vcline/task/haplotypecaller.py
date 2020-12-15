#!/usr/bin/env python

import re
from itertools import chain
from pathlib import Path

import luigi
from ftarc.task.base import ShellTask
from ftarc.task.picard import CreateSequenceDictionary
from ftarc.task.resource import (FetchDbsnpVCF, FetchMillsIndelVCF,
                                 FetchReferenceFASTA)
from ftarc.task.samtools import samtools_view_and_index
from luigi.util import requires

from .align import PrepareCRAMNormal
from .ref import FetchEvaluationIntervalList, FetchHapmapVCF
from .samtools import samtools_merge_and_index


@requires(FetchEvaluationIntervalList, FetchReferenceFASTA,
          CreateSequenceDictionary)
class SplitEvaluationIntervals(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        scatter_count = self.cf['n_worker']
        if scatter_count > 1:
            run_dir = Path(self.cf['qc_dir_path']).joinpath(
                'intervals/{0}.split_in_{1}'.format(
                    Path(self.input()[0].path).stem, scatter_count
                )
            )
            return [
                luigi.LocalTarget(
                    run_dir.joinpath(f'{i:04d}-scattered.interval_list')
                ) for i in range(scatter_count)
            ]
        else:
            return [luigi.LocalTarget(self.input()[0].path)]

    def run(self):
        input_interval = Path(self.input()[0].path)
        run_id = input_interval.stem
        output_interval_paths = [o.path for o in self.output()]
        scatter_count = len(output_interval_paths)
        self.print_log(f'Split an interval list in {scatter_count}:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        fa_path = self.input()[1][0].path
        run_dir = Path(output_interval_paths[0]).parent
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=run_dir, remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} SplitIntervals'
                + f' --reference {fa_path}'
                + f' --intervals {input_interval}'
                + f' --scatter-count {scatter_count}'
                + f' --output {run_dir}'
            ),
            input_files_or_dirs=[input_interval, fa_path],
            output_files_or_dirs=[*output_interval_paths, run_dir]
        )


@requires(PrepareCRAMNormal, FetchReferenceFASTA, FetchDbsnpVCF,
          SplitEvaluationIntervals, CreateSequenceDictionary)
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
        interval_paths = [i.path for i in self.input()[3]]
        skip_interval_split = (len(interval_paths) == 1)
        fa_path = self.input()[1][0].path
        input_cram_path = self.input()[0][0].path
        dbsnp_vcf_path = self.input()[2][0].path
        output_path_prefix = '.'.join(str(output_gvcf).split('.')[:-3])
        if skip_interval_split:
            tmp_prefixes = [output_path_prefix]
        else:
            tmp_prefixes = [
                '{0}.{1}'.format(output_path_prefix, Path(p).stem)
                for p in interval_paths
            ]
        input_targets = yield [
            HaplotypeCaller(
                input_cram_path=input_cram_path, fa_path=fa_path,
                dbsnp_vcf_path=dbsnp_vcf_path, evaluation_interval_path=p,
                output_path_prefix=s, cf=self.cf
            ) for p, s in zip(interval_paths, tmp_prefixes)
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
        if skip_interval_split:
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
            tmp_file_paths = list(
                chain.from_iterable([
                    [o.path for o in t] for t in input_targets
                ])
            )
            self.run_shell(
                args=('rm -f' + ''.join([f' {p}' for p in tmp_file_paths])),
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
        run_dir = output_gvcf.parent
        self.setup_shell(
            run_id='.'.join(output_gvcf.name.split('.')[:-3]),
            log_dir_path=self.cf['log_dir_path'], commands=gatk, cwd=run_dir,
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
            output_files_or_dirs=[*output_file_paths, run_dir]
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
