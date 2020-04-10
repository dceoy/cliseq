#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareNormalCRAM, PrepareTumorCRAM
from .base import ShellTask
from .haplotypecaller import GenotypeGVCF
from .ref import (CreateSequenceDictionary, FetchCnvBlackList,
                  FetchEvaluationIntervalList, FetchReferenceFASTA)


@requires(GenotypeGVCF)
class CreateGermlineSnpIntervalList(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
                    Path(self.input()[0].path).stem + '.interval_list'
                )
            )
        )

    def run(self):
        input_vcf_path = self.input()[0].path
        run_id = Path(input_vcf_path).stem
        self.print_log(f'Create a germline SNP interval_list:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        output_interval_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_cnv_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} VcfToIntervalList'
                + f' --INPUT {input_vcf_path}'
                + ' --INCLUDE_FILTERED true'
                + f' --OUTPUT {output_interval_path}'
            ),
            input_files_or_dirs=input_vcf_path,
            output_files_or_dirs=output_interval_path
        )


@requires(CreateGermlineSnpIntervalList, FetchReferenceFASTA)
class CollectAllelicCounts(ShellTask):
    cram_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
                    Path(self.cram_path).stem + '.allelic_counts.tsv'
                )
            )
        )

    def run(self):
        run_id = Path(self.cram_path).stem
        self.print_log(f'Collects allele counts:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        common_sites_interval_path = self.input()[0].path
        fa_path = self.input()[1].path
        allelic_counts_tsv_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_cnv_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CollectAllelicCounts'
                + f' --intervals {common_sites_interval_path}'
                + f' --input {self.cram_path}'
                + f' --reference {fa_path}'
                + ' --minimum-base-quality 20'
                + f' --output {allelic_counts_tsv_path}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                self.cram_path, common_sites_interval_path, fa_path
            ],
            output_files_or_dirs=allelic_counts_tsv_path
        )


@requires(PrepareTumorCRAM, PrepareNormalCRAM)
class CollectMatchedAllelicCounts(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            CollectAllelicCounts(
                cram_path=i[0].path, fq_list=self.fq_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            ).output() for i in self.input()
        ]


@requires(FetchEvaluationIntervalList, FetchCnvBlackList,
          FetchReferenceFASTA, CreateSequenceDictionary)
class PreprocessIntervals(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
                    Path(self.input()[0].path).stem
                    + '.preprocessed.interval_list'
                )
            )
        )

    def run(self):
        evaluation_interval_path = self.input()[0].path
        run_id = Path(evaluation_interval_path).stem
        self.print_log(f'Prepares bins for coverage collection:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        cnv_black_list_path = self.input()[1].path
        fa_path = self.input()[2].path
        seq_dict_path = self.input()[3].path
        preprocessed_interval_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_cnv_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} PreprocessIntervals'
                + f' --intervals {evaluation_interval_path}'
                + f' --exclude-intervals {cnv_black_list_path}'
                + f' --sequence-dictionary {seq_dict_path}'
                + f' --reference {fa_path}'
                + (
                    ' --bin-length 0 --padding 250' if self.cf['exome']
                    else ' --bin-length 1000 --padding 0'
                )
                + ' --interval-merging-rule OVERLAPPING_ONLY'
                + f' --output {preprocessed_interval_path}'
            ),
            input_files_or_dirs=[
                evaluation_interval_path, cnv_black_list_path,
                seq_dict_path, fa_path
            ],
            output_files_or_dirs=preprocessed_interval_path
        )


@requires(PrepareTumorCRAM, PreprocessIntervals, FetchReferenceFASTA)
class CollectReadCountsTumor(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
                    Path(self.input()[0][0].path).stem + '.counts.hdf5'
                )
            )
        )

    def run(self):
        input_cram_path = self.input()[0][0].path
        run_id = Path(input_cram_path).stem
        self.print_log(f'Collects read counts:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        preprocessed_interval_path = self.input()[1].path
        fa_path = self.input()[2].path
        counts_hdf5_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_cnv_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CollectReadCounts'
                + f' --intervals {preprocessed_interval_path}'
                + f' --input {input_cram_path}'
                + f' --reference {fa_path}'
                + ' --format HDF5'
                + ' --interval-merging-rule OVERLAPPING_ONLY'
                + f' --output {counts_hdf5_path}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                input_cram_path, preprocessed_interval_path, fa_path
            ],
            output_files_or_dirs=counts_hdf5_path
        )


@requires(CollectReadCountsTumor)
class DenoiseReadCountsTumor(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
                        Path(Path(self.input().path).stem).stem + f'.{s}.tsv'
                    )
                )
            ) for s in ['denoised_cr', 'standardized_cr']
        ]

    def run(self):
        input_counts_path = self.input().path
        run_id = Path(input_counts_path).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        denoised_cr_tsv_path = self.output()[0].path
        standardized_cr_tsv_path = self.output()[1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_cnv_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} DenoiseReadCounts'
                + f' --input {input_counts_path}'
                + f' --standardized-copy-ratios {standardized_cr_tsv_path}'
                + f' --denoised-copy-ratios {denoised_cr_tsv_path}'
            ),
            input_files_or_dirs=[input_counts_path],
            output_files_or_dirs=[
                standardized_cr_tsv_path, denoised_cr_tsv_path
            ]
        )


@requires(DenoiseReadCountsTumor, CollectMatchedAllelicCounts)
class ModelSegmentsTumor(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
                        create_matched_id(*[i.path for i in self.input()[1]])
                        + f'.{s}'
                    )
                )
            ) for s in [
                'cr.seg', 'hets.tsv', 'hets.normal.tsv', 'cr.igv.seg',
                'af.igv.seg', 'modelBegin.seg', 'modelBegin.cr.param',
                'modelBegin.af.param', 'modelFinal.seg', 'modelFinal.cr.param',
                'modelFinal.af.param'
            ]
        ]

    def run(self):
        output_dir_paths = [o.path for o in self.output()]
        run_id = Path(Path(output_dir_paths[0]).stem).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        denoised_cr_tsv_path = self.input()[0][0].path
        allelic_counts_tsv_path = self.input()[1][0].path
        normal_allelic_counts_tsv_path = self.input()[1][1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_cnv_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} ModelSegments'
                + f' --denoised-copy-ratios {denoised_cr_tsv_path}'
                + f' --allelic-counts {allelic_counts_tsv_path}'
                + f' --normal-allelic-counts {normal_allelic_counts_tsv_path}'
                + ' --minimum-total-allele-count-case 0'
                + ' --minimum-total-allele-count-normal 30'
                + ' --genotyping-homozygous-log-ratio-threshold -10.0'
                + ' --genotyping-base-error-rate 0.05'
                + ' --maximum-number-of-segments-per-chromosome 1000'
                + ' --kernel-variance-copy-ratio 0.0'
                + ' --kernel-variance-allele-fraction 0.025'
                + ' --kernel-scaling-allele-fraction 1.0'
                + ' --kernel-approximation-dimension 100'
                + ''.join([
                    f' --window-size {v}' for v in [8, 16, 32, 64, 128, 256]
                ])
                + ' --number-of-changepoints-penalty-factor 1.0'
                + ' --minor-allele-fraction-prior-alpha 25.0'
                + ' --number-of-samples-copy-ratio 100'
                + ' --number-of-burn-in-samples-copy-ratio 50'
                + ' --number-of-samples-allele-fraction 100'
                + ' --number-of-burn-in-samples-allele-fraction 50'
                + ' --smoothing-credible-interval-threshold-copy-ratio 2.0'
                + ' --smoothing-credible-interval-threshold-allele-fraction'
                + ' 2.0'
                + ' --maximum-number-of-smoothing-iterations 10'
                + ' --number-of-smoothing-iterations-per-fit 0'
                + ' --output {}'.format(
                    self.cf['somatic_cnv_gatk_dir_path']
                )
                + ' --output-prefix {run_id}'
            ),
            input_files_or_dirs=[
                denoised_cr_tsv_path, allelic_counts_tsv_path,
                normal_allelic_counts_tsv_path
            ],
            output_files_or_dirs=output_dir_paths
        )


@requires(ModelSegmentsTumor)
class CallCopyRatioSegmentsTumor(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
                        Path(self.input()[0].path).stem + '.called.seg'
                    )
                )
            )
        ]

    def run(self):
        cr_seg_path = self.input()[0].path
        run_id = cr_seg_path
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        output_seg_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['somatic_cnv_gatk_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CallCopyRatioSegments'
                + f' --input {cr_seg_path}'
                + ' --neutral-segment-copy-ratio-lower-bound 0.9'
                + ' --neutral-segment-copy-ratio-upper-bound 1.1'
                + ' --outlier-neutral-segment-copy-ratio-z-score-threshold 2.0'
                + ' --calling-copy-ratio-z-score-threshold 2.0'
                + f' --output {output_seg_path}'
            ),
            input_files_or_dirs=cr_seg_path,
            output_files_or_dirs=output_seg_path
        )


if __name__ == '__main__':
    luigi.run()
