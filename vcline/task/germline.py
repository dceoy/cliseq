#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .align import PrepareCRAMs
from .base import ShellTask
from .ref import (CreateEvaluationIntervalList, CreateFASTAIndex,
                  FetchDbsnpVCF, FetchReferenceFASTA,
                  PrepareGermlineResourceVCFs)


@requires(CreateEvaluationIntervalList, FetchReferenceFASTA)
class SplitEvaluationIntervals(ShellTask):
    cf = luigi.DictParameter()
    priority = 40

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['call_dir_path']).joinpath(
                        f'{i:04d}-scattered.interval_list'
                    )
                )
            ) for i in range(self.cf['n_cpu_per_worker'])
        ]

    def run(self):
        interval_path = self.input()[0].path
        run_id = Path(interval_path).stem
        print_log(f'Split an evaluation interval list:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        fa_path = self.input()[1].path
        scatter_count = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['call_dir_path']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} SplitIntervals'
                + f' --reference {fa_path}'
                + f' --intervals {interval_path}'
                + f' --scatter-count {scatter_count}'
                + ' --output {}'.format(self.cf['call_dir_path'])
            ),
            input_files=[interval_path, fa_path],
            output_files=[o.path for o in self.output()]
        )


class PrepareEvaluationIntervals(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    evaluation_interval_path = luigi.Parameter(default='')
    cf = luigi.DictParameter()
    priority = 40

    def requires(self):
        return (
            SplitEvaluationIntervals(
                ref_fa_paths=self.ref_fa_paths,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            ) if self.cf['split_intervals'] else
            CreateEvaluationIntervalList(
                ref_fa_paths=self.ref_fa_paths,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        )

    def output(self):
        return (
            self.input() if isinstance(self.input(), list) else [self.input()]
        )


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          FetchDbsnpVCF, PrepareEvaluationIntervals)
class CallVariantsWithHaplotypeCaller(ShellTask):
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['call_dir_path']).joinpath(
                        Path(self.input()[0][1][0].path).stem
                        + f'.HaplotypeCaller.{s}'
                    )
                )
            ) for s in ['raw.g.vcf.gz', 'cram', 'cram.crai']
        ]

    def run(self):
        gvcf_path = self.output()[0].path
        run_id = '.'.join(Path(gvcf_path).name.split('.')[:-4])
        print_log(f'Call germline variants with HaplotypeCaller:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        samtools = self.cf['samtools']
        save_memory = str(self.cf['memory_mb_per_worker'] < 8 * 1024).lower()
        n_cpu = self.cf['n_cpu_per_worker']
        memory_per_thread = self.cf['samtools_memory_per_thread']
        input_cram_path = self.input()[0][1][0].path
        fa_path = self.input()[1].path
        fai_path = self.input()[2].path
        dbsnp_vcf_path = self.input()[3][0].path
        evaluation_interval_paths = [i.path for i in self.input()[4]]
        output_cram_path = self.output()[1].path
        tmp_bam_paths = [
            re.sub(
                r'(\.cram)$', '.{}.bam'.format(Path(i).stem), output_cram_path
            ) for i in evaluation_interval_paths
        ]
        tmp_gvcf_paths = (
            [
                re.sub(
                    r'\.g\.vcf\.gz$', '.{}.g.vcf.gz'.format(Path(i).stem),
                    gvcf_path
                ) for i in evaluation_interval_paths
            ] if len(evaluation_interval_paths) > 1
            else [gvcf_path]
        )
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, samtools], cwd=self.cf['call_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {gatk}{gatk_opts} HaplotypeCaller'
                    + f' --reference {fa_path}'
                    + f' --input {input_cram_path}'
                    + f' --dbsnp {dbsnp_vcf_path}'
                    + f' --intervals {i}'
                    + f' --output {g}'
                    + f' --bam-output {b}'
                    + ' --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP'
                    + f' --native-pair-hmm-threads {n_cpu}'
                    + f' --disable-bam-index-caching {save_memory}'
                    + ' --emit-ref-confidence GVCF'
                    + ''.join(
                        [
                            f' --annotation-group {g}' for g in [
                                'StandardAnnotation', 'AS_StandardAnnotation',
                                'StandardHCAnnotation'
                            ]
                        ] + [
                            f' --gvcf-gq-bands {i}' for i in range(10, 100, 10)
                        ]
                    )
                    + ' --create-output-bam-index false'
                ) for i, g, b in zip(
                    evaluation_interval_paths, tmp_gvcf_paths, tmp_bam_paths
                )
            ],
            input_files=[
                input_cram_path, fa_path, dbsnp_vcf_path,
                *evaluation_interval_paths
            ],
            output_files=[*tmp_gvcf_paths, *tmp_bam_paths],
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
        if len(tmp_gvcf_paths) > 1:
            self.run_shell(
                args=[
                    (
                        f'set -e && {gatk}{gatk_opts} CombineGVCFs'
                        + f' --reference {fa_path}'
                        + ''.join([f' --variant {g}' for g in tmp_gvcf_paths])
                        + f' --output {gvcf_path}'
                    ),
                    (f'set -e && rm -f ' + ' '.join(tmp_gvcf_paths))
                ],
                input_files=[*tmp_gvcf_paths, fa_path], output_files=gvcf_path
            )


@requires(CallVariantsWithHaplotypeCaller, FetchReferenceFASTA,
          FetchDbsnpVCF, PrepareEvaluationIntervals)
class GenotypeGVCF(ShellTask):
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.g\.vcf\.gz$', '.vcf.gz', self.input()[0][0].path)
        )

    def run(self):
        vcf_path = self.output().path
        run_id = '.'.join(Path(vcf_path).name.split('.')[:-3])
        print_log(f'Genotype a HaplotypeCaller GVCF:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['memory_mb_per_worker'] < 8 * 1024).lower()
        gvcf_path = self.input()[0].path
        fa_path = self.input()[1].path
        dbsnp_vcf_path = self.input()[2][0].path
        evaluation_interval_path = self.input()[3].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['call_dir_path']
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
            input_files=[
                gvcf_path, fa_path, dbsnp_vcf_path, evaluation_interval_path
            ],
            output_files=vcf_path
        )


@requires(GenotypeGVCF, FetchReferenceFASTA, PrepareGermlineResourceVCFs)
class ApplyVQSR(ShellTask):
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.raw\.vcf\.gz$', '.vcf.gz', self.input()[0].path)
        )

    def run(self):
        filtered_vcf_path = self.output().path
        run_id = '.'.join(Path(filtered_vcf_path).name.split('.')[:-2])
        print_log(f'Apply Variant Quality Score Recalibration:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['memory_mb_per_worker'] < 8 * 1024).lower()
        raw_vcf_path = self.input()[0].path
        fa_path = self.input()[1].path
        resource_params = {
            'hapmap': 'hapmap,known=false,training=true,truth=true,prior=15.0',
            'omni': 'omni,known=false,training=true,truth=false,prior=12.0',
            'snp_1000G':
            '1000G,known=false,training=true,truth=false,prior=10.0',
            'dbsnp': 'dbsnp,known=true,training=false,truth=false,prior=2.0'
        }
        resource_vcf_paths = {k: v[0].path for k, v in self.input()[2].items()}
        recal_path = re.sub(r'\.vcf\.gz$', '.recal', filtered_vcf_path)
        tranches_path = re.sub(r'\.vcf\.gz$', '.tranches', filtered_vcf_path)
        plot_r_path = re.sub(r'\.vcf\.gz$', '.plot.R', filtered_vcf_path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['call_dir_path']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} VariantRecalibrator'
                + f' --reference {fa_path}'
                + f' --variant {raw_vcf_path}'
                + ''.join([
                    f' --resource:{resource_params[k]} {v}'
                    for k, v in resource_vcf_paths.items()
                ])
                + f' --output {recal_path}'
                + f' --tranches-file {tranches_path}'
                + f' --rscript-file {plot_r_path}'
                + ''.join([
                    f' --use-annotation {a}' for a in
                    ['QD', 'MQ', 'MQRankSum', 'ReadPosRankSum', 'FS', 'SOR']
                ])
                + ' --mode SNP'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files=[raw_vcf_path, fa_path, *resource_vcf_paths.values()],
            output_files=[recal_path, tranches_path, plot_r_path]
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} ApplyVQSR'
                + f' --reference {fa_path}'
                + f' --variant {raw_vcf_path}'
                + f' --tranches-file {tranches_path}'
                + f' --recal-file {recal_path}'
                + f' --output {filtered_vcf_path}'
                + ' --truth-sensitivity-filter-level 99.0'
                + ' --mode SNP'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files=[raw_vcf_path, fa_path, recal_path, tranches_path],
            output_files=filtered_vcf_path
        )


if __name__ == '__main__':
    luigi.run()
