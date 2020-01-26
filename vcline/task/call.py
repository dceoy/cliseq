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


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          FetchEvaluationIntervalList)
class CallSomaticVariantsWithMutect2(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['call_dir_path']).joinpath(
                    create_matched_id(*[i[0].path for i in self.input()[0]])
                    + '.Mutect2.somtic.vcf'
                )
            )
        )

    def run(self):
        run_id = '.'.join(Path(self.output().path).stem.split('.')[:-2])
        print_log(f'Call somatic short variants with Mutect2:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['memory_mb_per_worker'] < 8 * 1024).lower()
        cram_paths = [i[0].path for i in self.input()[0]]
        fa_path = self.input()[1].path
        evaluation_interval_path = self.input()[3].path
        filtered_vcf_path = self.output().path
        unfiltered_vcf_path = re.sub(
            r'(\.vcf)$', '.unfiltered\\1', filtered_vcf_path
        )
        normal_name = self.sample_names[1]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['call_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} Mutect2'
                + f' --reference {fa_path}'
                + f' --input {cram_paths[0]}'
                + f' --input {cram_paths[1]}'
                + f' --output {unfiltered_vcf_path}'
                + f' --normal-sample {normal_name}'
                + f' --intervals {evaluation_interval_path}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files=[*cram_paths, fa_path],
            output_files=unfiltered_vcf_path
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} FilterMutectCalls'
                + f' --reference {fa_path}'
                + f' --variant {unfiltered_vcf_path}'
                + f' --output {filtered_vcf_path}'
            ),
            input_files=[unfiltered_vcf_path, fa_path],
            output_files=filtered_vcf_path
        )


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          FetchDbsnpVCF, FetchEvaluationIntervalList)
class CallGermlineVariantsWithHaplotypeCaller(ShellTask):
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['call_dir_path']).joinpath(
                    Path(self.input()[0][1][0].path).stem
                    + '.HaplotypeCaller.gernline.vcf'
                )
            )
        )

    def run(self):
        run_id = '.'.join(Path(self.output().path).stem.split('.')[:-2])
        print_log(
            f'Call germline SNPs and indels with HaplotypeCaller:\t{run_id}'
        )
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['memory_mb_per_worker'] < 8 * 1024).lower()
        cram_path = self.input()[0][1][0].path
        vcf_gz_path = self.output().path
        gvcf_gz_path = re.sub(r'\.vcf$', '.g.vcf.gz', vcf_gz_path)
        fa_path = self.input()[1].path
        dbsnp_vcf_path = self.input()[3][0].path
        evaluation_interval_path = self.input()[4].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['call_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} HaplotypeCaller'
                + f' --reference {fa_path}'
                + f' --input {cram_path}'
                + f' --dbsnp {dbsnp_vcf_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {gvcf_gz_path}'
                + f' --disable-bam-index-caching {save_memory}'
                + ' --emit-ref-confidence GVCF'
            ),
            input_files=[
                cram_path, fa_path, dbsnp_vcf_path, evaluation_interval_path
            ],
            output_files=gvcf_gz_path
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} GenotypeGVCFs'
                + f' --reference {fa_path}'
                + f' --variant {gvcf_gz_path}'
                + f' --dbsnp {dbsnp_vcf_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {vcf_gz_path}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files=[
                gvcf_gz_path, fa_path, dbsnp_vcf_path, evaluation_interval_path
            ],
            output_files=vcf_gz_path
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
            CallSomaticVariantsWithMutect2(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            ),
            CallGermlineVariantsWithHaplotypeCaller(
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
