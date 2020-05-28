#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .base import ShellTask
from .haplotypecaller import GenotypeHaplotypeCallerGVCF
from .ref import (CreateSequenceDictionary, FetchReferenceFASTA,
                  PreprocessIntervals)


@requires(GenotypeHaplotypeCallerGVCF)
class CreateGermlineSnpIntervalList(ShellTask):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['qc_dir_path']).joinpath('cnv').joinpath(
                Path(self.input()[0].path).stem + '.interval_list'
            )
        )

    def run(self):
        input_vcf = Path(self.input()[0].path)
        run_id = input_vcf.stem
        self.print_log(f'Create a germline SNP interval_list:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        output_interval = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=output_interval.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} VcfToIntervalList'
                + f' --INPUT {input_vcf}'
                + ' --INCLUDE_FILTERED true'
                + f' --OUTPUT {output_interval}'
            ),
            input_files_or_dirs=input_vcf, output_files_or_dirs=output_interval
        )


class CollectAllelicCounts(ShellTask):
    cram_path = luigi.Parameter()
    common_sites_interval_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['qc_dir_path']).joinpath('cnv').joinpath(
                Path(self.cram_path).stem + '.allelic_counts.tsv'
            )
        )

    def run(self):
        run_id = Path(self.cram_path).stem
        self.print_log(f'Collects allele counts:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        allelic_counts_tsv = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=allelic_counts_tsv.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CollectAllelicCounts'
                + f' --intervals {self.common_sites_interval_path}'
                + f' --input {self.cram_path}'
                + f' --reference {self.fa_path}'
                + f' --output {allelic_counts_tsv}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                self.cram_path, self.common_sites_interval_path, self.fa_path
            ],
            output_files_or_dirs=allelic_counts_tsv
        )


@requires(PrepareCRAMTumor, CreateGermlineSnpIntervalList, FetchReferenceFASTA,
          CreateSequenceDictionary)
class CollectAllelicCountsTumor(luigi.Task):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['qc_dir_path']).joinpath('cnv').joinpath(
                Path(self.input()[0][0].path).stem + '.allelic_counts.tsv'
            )
        )

    def run(self):
        yield CollectAllelicCounts(
            cram_path=self.input()[0][0].path,
            common_sites_interval_path=self.input()[1].path,
            fa_path=self.input()[2][0].path, cf=self.cf
        )


@requires(PrepareCRAMNormal, CreateGermlineSnpIntervalList,
          FetchReferenceFASTA, CreateSequenceDictionary)
class CollectAllelicCountsNormal(luigi.Task):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['qc_dir_path']).joinpath('cnv').joinpath(
                Path(self.input()[0][0].path).stem + '.allelic_counts.tsv'
            )
        )

    def run(self):
        yield CollectAllelicCounts(
            cram_path=self.input()[0][0].path,
            common_sites_interval_path=self.input()[1].path,
            fa_path=self.input()[2][0].path, cf=self.cf
        )


class CollectReadCounts(ShellTask):
    cram_path = luigi.Parameter()
    preprocessed_interval_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['qc_dir_path']).joinpath('cnv').joinpath(
                Path(self.cram_path).stem + '.counts.hdf5'
            )
        )

    def run(self):
        run_id = Path(self.cram_path).stem
        self.print_log(f'Collects read counts:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        save_memory = str(self.cf['save_memory']).lower()
        counts_hdf5 = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=counts_hdf5.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CollectReadCounts'
                + f' --intervals {self.preprocessed_interval_path}'
                + f' --input {self.cram_path}'
                + f' --reference {self.fa_path}'
                + ' --format HDF5'
                + ' --interval-merging-rule OVERLAPPING_ONLY'
                + f' --output {counts_hdf5}'
                + f' --disable-bam-index-caching {save_memory}'
            ),
            input_files_or_dirs=[
                self.cram_path, self.preprocessed_interval_path, self.fa_path
            ],
            output_files_or_dirs=counts_hdf5
        )


@requires(CollectReadCounts)
class DenoiseReadCounts(ShellTask):
    seq_dict_path = luigi.Parameter()
    cf = luigi.DictParameter()
    create_plots = luigi.BoolParameter(default=True)
    priority = 30

    def output(self):
        output_path_prefix = re.sub(r'\.counts\.hdf5$', '', self.input().path)
        return [
            luigi.LocalTarget(f'{output_path_prefix}.{s}')
            for s in ['denoised_cr.tsv', 'standardized_cr.tsv']
        ]

    def run(self):
        counts_hdf5_path = self.input().path
        run_id = Path(Path(counts_hdf5_path).stem).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        r = self.cf['R']
        denoised_cr_tsv_path = self.output()[0].path
        standardized_cr_tsv_path = self.output()[1].path
        run_dir = Path(denoised_cr_tsv_path).parent
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, r], cwd=run_dir,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} DenoiseReadCounts'
                + f' --input {counts_hdf5_path}'
                + f' --standardized-copy-ratios {standardized_cr_tsv_path}'
                + f' --denoised-copy-ratios {denoised_cr_tsv_path}'
            ),
            input_files_or_dirs=counts_hdf5_path,
            output_files_or_dirs=[
                standardized_cr_tsv_path, denoised_cr_tsv_path
            ]
        )
        if self.create_plots:
            plots_dir = run_dir.joinpath(f'{run_id}.plots')
            self.run_shell(
                args=(
                    f'set -e && {gatk}{gatk_opts} PlotDenoisedCopyRatios'
                    + f' --standardized-copy-ratios {standardized_cr_tsv_path}'
                    + f' --denoised-copy-ratios {denoised_cr_tsv_path}'
                    + f' --sequence-dictionary {self.seq_dict_path}'
                    + f' --output {plots_dir}'
                    + f' --output-prefix {run_id}'

                ),
                input_files_or_dirs=[
                    standardized_cr_tsv_path, denoised_cr_tsv_path,
                    self.seq_dict_path
                ],
                output_files_or_dirs=[
                    plots_dir.joinpath(f'{run_id}.denoised.png'), plots_dir
                ]
            )


@requires(DenoiseReadCounts)
class ModelSegments(ShellTask):
    seq_dict_path = luigi.Parameter()
    case_allelic_counts_tsv_path = luigi.Parameter(default='')
    normal_allelic_counts_tsv_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter()
    cf = luigi.DictParameter()
    min_total_allele_counts = luigi.DictParameter(
        default={'case': 0, 'normal': 30}
    )
    create_plots = luigi.BoolParameter(default=True)
    priority = 30

    def output(self):
        run_dir = Path(self.dest_dir_path)
        return [
            luigi.LocalTarget(run_dir.joinpath(f'{run_dir.name}.{s}'))
            for s in ['cr.seg', 'hets.tsv', 'modelFinal.seg']
        ]

    def run(self):
        output_file_paths = [o.path for o in self.output()]
        run_id = Path(Path(output_file_paths[0]).stem).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        r = self.cf['R']
        denoised_cr_tsv_path = self.input()[0].path
        run_dir = Path(output_file_paths[0]).parent
        if self.case_allelic_counts_tsv_path:
            allelic_count_args = (
                f' --allelic-counts {self.case_allelic_counts_tsv_path}'
                + ' --normal-allelic-counts'
                + f' {self.normal_allelic_counts_tsv_path}'
                + ''.join([
                    f' --minimum-total-allele-count-{k} {v}'
                    for k, v in self.min_total_allele_counts.items()
                ])
            )
            input_file_paths = [
                denoised_cr_tsv_path, self.case_allelic_counts_tsv_path,
                self.normal_allelic_counts_tsv_path
            ]
        else:
            allelic_count_args = (
                f' --allelic-counts {self.normal_allelic_counts_tsv_path}'
                + ' --minimum-total-allele-count-case {}'.format(
                    self.min_total_allele_counts['normal']
                )
            )
            input_file_paths = [
                denoised_cr_tsv_path, self.normal_allelic_counts_tsv_path
            ]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, r], cwd=run_dir,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} ModelSegments'
                + f' --denoised-copy-ratios {denoised_cr_tsv_path}'
                + allelic_count_args
                + f' --output-prefix {run_id}'
                + f' --output {run_dir}'
            ),
            input_files_or_dirs=input_file_paths,
            output_files_or_dirs=[*output_file_paths, run_dir]
        )
        if self.create_plots:
            plots_dir = run_dir.joinpath(f'{run_id}.plots')
            het_allelic_counts_tsv_path = output_file_paths[1]
            modeled_segments_path = output_file_paths[2]
            self.run_shell(
                args=(
                    f'set -e && {gatk}{gatk_opts} PlotModeledSegments'
                    + f' --denoised-copy-ratios {denoised_cr_tsv_path}'
                    + f' --allelic-counts {het_allelic_counts_tsv_path}'
                    + f' --segments {modeled_segments_path}'
                    + f' --sequence-dictionary {self.seq_dict_path}'
                    + f' --output {plots_dir}'
                    + f' --output-prefix {run_id}'

                ),
                input_files_or_dirs=[
                    denoised_cr_tsv_path, het_allelic_counts_tsv_path,
                    modeled_segments_path, self.seq_dict_path
                ],
                output_files_or_dirs=[
                    plots_dir.joinpath(f'{run_id}.modeled.png'), plots_dir
                ]
            )


@requires(ModelSegments)
class CallCopyRatioSegments(ShellTask):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.seg$', '.called.seg', self.input()[0].path)
        )

    def run(self):
        cr_seg_path = self.input()[0].path
        run_id = Path(cr_seg_path).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        output_seg = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=output_seg.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} CallCopyRatioSegments'
                + f' --input {cr_seg_path}'
                + f' --output {output_seg}'
            ),
            input_files_or_dirs=cr_seg_path, output_files_or_dirs=output_seg
        )


@requires(PrepareCRAMTumor, PreprocessIntervals, FetchReferenceFASTA,
          CreateSequenceDictionary, CollectAllelicCountsTumor,
          CollectAllelicCountsNormal)
class CallCopyRatioSegmentsTumor(luigi.Task):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        allelic_counts_paths = [i.path for i in self.input()[4:6]]
        run_dir = Path(allelic_counts_paths[0]).parent
        return luigi.LocalTarget(
            run_dir.joinpath(
                Path(create_matched_id(*allelic_counts_paths)).stem
                + '.cr.called.seg'
            )
        )

    def run(self):
        yield CallCopyRatioSegments(
            cram_path=self.input()[0][0].path,
            preprocessed_interval_path=self.input()[1].path,
            fa_path=self.input()[2][0].path,
            seq_dict_path=self.input()[3].path,
            case_allelic_counts_tsv_path=self.input()[4].path,
            normal_allelic_counts_tsv_path=self.input()[5].path,
            dest_dir_path=str(Path(self.output().path).parent), cf=self.cf
        )


@requires(PrepareCRAMNormal, PreprocessIntervals, FetchReferenceFASTA,
          CreateSequenceDictionary, CollectAllelicCountsNormal)
class CallCopyRatioSegmentsNormal(luigi.Task):
    cf = luigi.DictParameter()
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            Path(self.input()[4].path).parent.joinpath(
                Path(Path(self.input()[4].path).stem).stem + '.cr.called.seg'
            )
        )

    def run(self):
        yield CallCopyRatioSegments(
            cram_path=self.input()[0][0].path,
            preprocessed_interval_path=self.input()[1].path,
            fa_path=self.input()[2][0].path,
            seq_dict_path=self.input()[3].path,
            normal_allelic_counts_tsv_path=self.input()[4].path,
            dest_dir_path=str(Path(self.output().path).parent), cf=self.cf
        )


@requires(CallCopyRatioSegmentsTumor, CallCopyRatioSegmentsNormal)
class CallCopyRatioSegmentsMatched(luigi.WrapperTask):
    priority = 30

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
