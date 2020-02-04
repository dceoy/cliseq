#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id, print_log
from .align import PrepareCRAMs
from .base import ShellTask
from .ref import (CreateFASTAIndex, FetchDbsnpVCF, FetchEvaluationIntervalList,
                  FetchReferenceFASTA)


@requires(FetchEvaluationIntervalList, FetchReferenceFASTA)
class SplitEvaluationIntervals(ShellTask):
    cf = luigi.DictParameter()
    priority = 40

    def output(self):
        return (
            [
                luigi.LocalTarget(
                    str(
                        Path(self.cf['call_dir_path']).joinpath(
                            f'{i:04d}-scattered.interval_list'
                        )
                    )
                ) for i in range(self.cf['n_cpu_per_worker'])
            ] if self.cf['n_cpu_per_worker'] > 1 else
            [luigi.LocalTarget(self.input()[0].path)]
        )

    def run(self):
        n_cpu = self.cf['n_cpu_per_worker']
        if n_cpu > 1:
            interval_path = self.input()[0].path
            run_id = Path(interval_path).stem
            print_log(f'Split an interval list:\t{run_id}')
            fa_path = self.input()[1].path
            gatk = self.cf['gatk']
            gatk_opts = ' --java-options "{}"'.format(
                self.cf['gatk_java_options']
            )
            self.setup_shell(
                run_id=run_id, log_dir_path=self.cf['log_dir_path'],
                commands=gatk, cwd=self.cf['call_dir_path']
            )
            self.run_shell(
                args=(
                    'set -e && '
                    + f'{gatk}{gatk_opts} SplitIntervals'
                    + f' --reference {fa_path}'
                    + f' --intervals {interval_path}'
                    + f' --scatter-count {n_cpu}'
                    + ' --output {}'.format(self.cf['call_dir_path'])
                ),
                input_files=[interval_path, fa_path],
                output_files=[o.path for o in self.output()]
            )


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          FetchDbsnpVCF, SplitEvaluationIntervals)
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
            ) for s in ['g.vcf.gz', 'cram', 'cram.crai']
        ]

    def run(self):
        gvcf_path = self.output()[0].path
        run_id = '.'.join(Path(gvcf_path).stem.split('.')[:-4])
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
            output_files=[*tmp_gvcf_paths, *tmp_bam_paths], asynchronous=True
        )
        self.run_shell(
            args=[
                (
                    (
                        'set -eo pipefail && '
                        + f'{samtools} merge -@ {n_cpu} -h - '
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
          FetchDbsnpVCF, FetchEvaluationIntervalList)
class GenotypeHaplotypeCallerGVCF(ShellTask):
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.g\.vcf\.gz$', '.vcf', self.input()[0][0].path)
        )

    def run(self):
        vcf_path = self.output().path
        run_id = '.'.join(vcf_path.split('.')[:-2])
        print_log(f'Genotype a HaplotypeCaller GVCF:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['memory_mb_per_worker'] < 8 * 1024).lower()
        gvcf_path = self.input()[0].path
        fa_path = self.input()[1].path
        dbsnp_vcf_path = self.input()[2][0].path
        evaluation_interval_path = self.input()[3].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['call_dir_path']
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


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          SplitEvaluationIntervals)
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
            ) for s in ['raw.vcf.gz', 'raw.vcf.gz.stats', 'cram', 'cram.crai']
        ]

    def run(self):
        raw_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(raw_vcf_path).stem.split('.')[:-4])
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
        tmp_bam_paths = [
            re.sub(
                r'(\.cram)$', '.{}.bam'.format(Path(i).stem), output_cram_path
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
                    + f' --normal-sample {normal_name}'
                    + f' --disable-bam-index-caching {save_memory}'
                    + ' --max-mnp-distance 0'
                    + ' --create-output-bam-index false'
                ) for i, v, b in zip(
                    evaluation_interval_paths, tmp_vcf_paths, tmp_bam_paths
                )
            ],
            input_files=[
                *input_cram_paths, fa_path, *evaluation_interval_paths
            ],
            output_files=[*tmp_vcf_paths, *tmp_bam_paths, *tmp_stats_paths],
            asynchronous=True
        )
        self.run_shell(
            args=[
                (
                    (
                        'set -eo pipefail && '
                        + f'{samtools} merge -@ {n_cpu} -h - '
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
          FetchEvaluationIntervalList)
class FilterMutect2Calls(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.raw\.vcf\.gz$', '.vcf', self.input()[0][0].path)
        )

    def run(self):
        filtered_vcf_path = self.output().path
        run_id = '.'.join(Path(filtered_vcf_path).stem.split('.')[:-2])
        print_log(f'Filter somatic variants called by Mutect2:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        raw_vcf_path = self.input()[0].path
        fa_path = self.input()[1].path
        evaluation_interval_path = self.input()[2].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['call_dir_path']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} FilterMutectCalls'
                + f' --reference {fa_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --variant {raw_vcf_path}'
                + f' --output {filtered_vcf_path}'
            ),
            input_files=[raw_vcf_path, fa_path, evaluation_interval_path],
            output_files=filtered_vcf_path
        )


class CallVariants(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    known_indel_vcf_paths = luigi.ListParameter()
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        return [
            GenotypeHaplotypeCallerGVCF(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
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
