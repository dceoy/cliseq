#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id, print_log
from .align import PrepareCRAMs
from .base import ShellTask
from .ref import (CreateEvaluationIntervalList, CreateFASTAIndex,
                  FetchDbsnpVCF, FetchEvaluationIntervalList,
                  FetchReferenceFASTA, FetchResourceVCF)


class SplitEvaluationIntervals(ShellTask):
    ref_fa_paths = luigi.ListParameter()
    evaluation_interval_path = luigi.Parameter(default='')
    cf = luigi.DictParameter()
    priority = 40

    def requires(self):
        return [
            FetchReferenceFASTA(
                ref_fa_paths=self.ref_fa_paths, cf=self.cf
            ),
            (
                FetchEvaluationIntervalList(
                    evaluation_interval_path=self.evaluation_interval_path,
                    cf=self.cf
                ) if self.evaluation_interval_path
                else CreateEvaluationIntervalList(
                    ref_fa_paths=self.ref_fa_paths, cf=self.cf
                )
            )
        ]

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
        fa_path = self.input()[0].path
        run_id = Path(fa_path).stem
        print_log(f'Split an evaluation interval list:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        interval_path = self.input()[1].path
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
        if self.cf['split_intervals']:
            return SplitEvaluationIntervals(
                ref_fa_paths=self.ref_fa_paths,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.evaluation_interval_path:
            return FetchEvaluationIntervalList(
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        else:
            return CreateEvaluationIntervalList(
                ref_fa_paths=self.ref_fa_paths, cf=self.cf
            )

    def output(self):
        return (
            self.input() if isinstance(self.input(), list) else [self.input()]
        )


class PrepareGermlineResourceVCFs(luigi.WrapperTask):
    dbsnp_vcf_path = luigi.Parameter()
    hapmap_vcf_path = luigi.Parameter()
    omni_vcf_path = luigi.Parameter()
    snp_1000g_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 40

    def requires(self):
        return {
            'dbsnp': FetchDbsnpVCF(
                dbsnp_vcf_path=self.dbsnp_vcf_path, cf=self.cf
            ),
            'hapmap': FetchResourceVCF(
                resource_vcf_path=self.hapmap_vcf_path, cf=self.cf
            ),
            'omni': FetchResourceVCF(
                resource_vcf_path=self.omni_vcf_path, cf=self.cf
            ),
            'snp_1000g': FetchResourceVCF(
                resource_vcf_path=self.snp_1000g_vcf_path, cf=self.cf
            )
        }

    def output(self):
        return self.input()


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
                ]) + f' --output {recal_path}'
                + f' --tranches-file {tranches_path}'
                + f' --rscript-file {plot_r_path}'
                + ''.join([
                    f' --use-annotation {a}' for a in
                    ['QD', 'MQ', 'MQRankSum', 'ReadPosRankSum', 'FS', 'SOR']
                ]) + ' --mode SNP'
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


@requires(PrepareCRAMs, FetchReferenceFASTA, PrepareEvaluationIntervals,
          FetchDbsnpVCF)
class CalculateContamination(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['call_dir_path']).joinpath(
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
        dbsnp_vcf_path = self.input()[3][0].path
        pileup_table_paths = [f'{p}.pileup.table' for p in input_cram_paths]
        segment_table_path = self.output()[1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['call_dir_path']
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {gatk}{gatk_opts} GetPileupSummaries'
                    + f' --reference {fa_path}'
                    + f' --input {c}'
                    + f' --variant {dbsnp_vcf_path}'
                    + f' --intervals {evaluation_interval_path}'
                    + f' --output {t}'
                ) for c, t in zip(input_cram_paths, pileup_table_paths)
            ],
            input_files=[
                *input_cram_paths, fa_path, evaluation_interval_path,
                dbsnp_vcf_path
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
          PrepareEvaluationIntervals)
class CallVariantsWithMutect2(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['call_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0]]
                        ) + f'.Mutect2.{s}'
                    )
                )
            ) for s in [
                'raw.vcf.gz', 'raw.vcf.gz.stats', 'cram', 'cram.crai',
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
        raw_stats_path = self.output()[1].path
        output_cram_path = self.output()[2].path
        ob_priors_path = self.output()[4].path
        tmp_bam_paths = [
            re.sub(
                r'\.cram$', '.{}.bam'.format(Path(i).stem), output_cram_path
            ) for i in evaluation_interval_paths
        ]
        f1r2_paths = [
            re.sub(
                r'\.cram$', '.{}.f1r2.tar.gz'.format(Path(i).stem),
                output_cram_path
            ) for i in evaluation_interval_paths
        ]
        tmp_vcf_paths = (
            [
                re.sub(
                    r'\.raw\.vcf\.gz$', '.{}.raw.vcf.gz'.format(Path(i).stem),
                    raw_vcf_path
                ) for i in evaluation_interval_paths
            ] if len(evaluation_interval_paths) > 1 else [raw_vcf_path]
        )
        tmp_stats_paths = [f'{v}.stats' for v in tmp_vcf_paths]
        normal_name = self.sample_names[1]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, samtools], cwd=self.cf['call_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {gatk}{gatk_opts} Mutect2'
                    + f' --reference {fa_path}'
                    + ''.join([f' --input {p}' for p in input_cram_paths])
                    + f' --intervals {i}'
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
                *input_cram_paths, fa_path, *evaluation_interval_paths
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
          PrepareEvaluationIntervals)
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
        # contamination_table_path = self.input()[3][0].path
        # segment_table_path = self.input()[3][1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['call_dir_path']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} FilterMutectCalls'
                + f' --reference {fa_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --variant {raw_vcf_path}'
                + f' --stats {raw_stats_path}'
                # + f' --contamination-table {contamination_table_path}'
                # + f' --tumor-segmentation {segment_table_path}'
                + f' --orientation-bias-artifact-priors {ob_priors_path}'
                + f' --output {filtered_vcf_path}'
            ),
            input_files=[
                raw_vcf_path, fa_path, evaluation_interval_path,
                raw_stats_path, ob_priors_path
                # contamination_table_path, segment_table_path,
            ],
            output_files=filtered_vcf_path
        )


class CallVariants(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    known_indel_vcf_paths = luigi.ListParameter()
    hapmap_vcf_path = luigi.Parameter()
    omni_vcf_path = luigi.Parameter()
    snp_1000g_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter(default='')
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
                snp_1000g_vcf_path=self.snp_1000g_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            ),
            FilterMutect2Calls(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        ]

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
