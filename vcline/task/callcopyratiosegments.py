#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from ftarc.task.resource import CreateSequenceDictionary, FetchReferenceFasta
from luigi.util import requires

from .core import VclineTask
from .cram import PrepareCramNormal, PrepareCramTumor
from .haplotypecaller import FilterVariantTranches
from .resource import FetchCnvBlackList, FetchEvaluationIntervalList


@requires(FetchEvaluationIntervalList, FetchCnvBlackList,
          FetchReferenceFasta, CreateSequenceDictionary)
class PreprocessIntervals(VclineTask):
    cf = luigi.DictParameter()
    param_dict = luigi.DictParameter(default=dict())
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        evaluation_interval = Path(self.input()[0].path)
        return luigi.LocalTarget(
            evaluation_interval.parent.joinpath(
                '{0}.preprocessed.w{1}s.interval_list'.format(
                    evaluation_interval.stem,
                    ('e' if self.cf['exome'] else 'g')
                )
            )
        )

    def run(self):
        preprocessed_interval = Path(self.output().path)
        run_id = preprocessed_interval.stem
        self.print_log(f'Prepare bins for coverage collection:\t{run_id}')
        gatk = self.cf['gatk']
        evaluation_interval = Path(self.input()[0].path)
        cnv_blacklist = Path(self.input()[1].path)
        fa = Path(self.input()[2][0].path)
        seq_dict = Path(self.input()[3].path)
        param_dict = (
            self.param_dict or (
                {'bin-length': 0, 'padding': 250} if self.cf['exome']
                else {'bin-length': 1000, 'padding': 0}
            )
        )
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=preprocessed_interval.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} PreprocessIntervals'
                + f' --intervals {evaluation_interval}'
                + f' --exclude-intervals {cnv_blacklist}'
                + f' --sequence-dictionary {seq_dict}'
                + f' --reference {fa}'
                + ''.join(f' --{k} {v}' for k, v in param_dict.items())
                + ' --interval-merging-rule OVERLAPPING_ONLY'
                + f' --output {preprocessed_interval}'
            ),
            input_files_or_dirs=[
                evaluation_interval, cnv_blacklist, seq_dict, fa
            ],
            output_files_or_dirs=preprocessed_interval
        )


@requires(FilterVariantTranches, FetchReferenceFasta, CreateSequenceDictionary)
class CreateCommonSnpIntervalList(VclineTask):
    cf = luigi.DictParameter()
    min_ab = luigi.FloatParameter(default=0.5)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['qc_dir_path']).joinpath('cnv').joinpath(
                Path(Path(self.input()[0][0].path).stem).stem
                + '.biallelic_snp.interval_list'
            )
        )

    def run(self):
        output_interval = Path(self.output().path)
        run_id = output_interval.stem
        self.print_log(f'Create a common SNP interval_list:\t{run_id}')
        input_vcf = Path(self.input()[0][0].path)
        fa = Path(self.input()[1][0].path)
        dest_dir = output_interval.parent
        biallelic_snp_vcf = dest_dir.joinpath(f'{output_interval.stem}.vcf.gz')
        gatk = self.cf['gatk']
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=dest_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} SelectVariants'
                + f' --variant {input_vcf}'
                + f' --reference {fa}'
                + f' --output {biallelic_snp_vcf}'
                + ' --select-type-to-include SNP'
                + ' --restrict-alleles-to BIALLELIC'
                + ' --exclude-filtered'
                + ' --lenient'
            ),
            input_files_or_dirs=[input_vcf, fa],
            output_files_or_dirs=[
                biallelic_snp_vcf, f'{biallelic_snp_vcf}.tbi'
            ]
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} VcfToIntervalList'
                + f' --INPUT {biallelic_snp_vcf}'
                + f' --OUTPUT {output_interval}'
            ),
            input_files_or_dirs=biallelic_snp_vcf,
            output_files_or_dirs=output_interval
        )


class CollectAllelicCounts(VclineTask):
    cram_path = luigi.Parameter()
    common_sites_interval_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
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
        allelic_counts_tsv = Path(self.output().path)
        cram = Path(self.cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        common_sites_interval = Path(self.common_sites_interval_path).resolve()
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=allelic_counts_tsv.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} CollectAllelicCounts'
                + f' --intervals {common_sites_interval}'
                + f' --input {cram}'
                + f' --reference {fa}'
                + f' --output {allelic_counts_tsv}'
                + ' --disable-bam-index-caching '
                + str(self.cf['save_memory']).lower()
            ),
            input_files_or_dirs=[cram, common_sites_interval, fa],
            output_files_or_dirs=allelic_counts_tsv
        )


@requires(PrepareCramTumor, CreateCommonSnpIntervalList, FetchReferenceFasta,
          CreateSequenceDictionary)
class CollectAllelicCountsTumor(luigi.Task):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
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
            fa_path=self.input()[2][0].path, cf=self.cf,
            n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            sh_config=self.sh_config
        )


@requires(PrepareCramNormal, CreateCommonSnpIntervalList, FetchReferenceFasta,
          CreateSequenceDictionary)
class CollectAllelicCountsNormal(luigi.Task):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
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
            fa_path=self.input()[2][0].path, cf=self.cf,
            n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            sh_config=self.sh_config
        )


class CollectReadCounts(VclineTask):
    cram_path = luigi.Parameter()
    preprocessed_interval_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
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
        cram = Path(self.cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        preprocessed_interval = Path(self.preprocessed_interval_path).resolve()
        gatk = self.cf['gatk']
        counts_hdf5 = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=counts_hdf5.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} CollectReadCounts'
                + f' --intervals {preprocessed_interval}'
                + f' --input {cram}'
                + f' --reference {fa}'
                + ' --format HDF5'
                + ' --interval-merging-rule OVERLAPPING_ONLY'
                + f' --output {counts_hdf5}'
                + ' --disable-bam-index-caching '
                + str(self.cf['save_memory']).lower()
            ),
            input_files_or_dirs=[cram, preprocessed_interval, fa],
            output_files_or_dirs=counts_hdf5
        )


@requires(CollectReadCounts)
class DenoiseReadCounts(VclineTask):
    seq_dict_path = luigi.Parameter()
    cf = luigi.DictParameter()
    create_plots = luigi.BoolParameter(default=True)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        output_path_prefix = re.sub(r'\.counts\.hdf5$', '', self.input().path)
        return [
            luigi.LocalTarget(f'{output_path_prefix}.{s}')
            for s in ['denoised_cr.tsv', 'standardized_cr.tsv']
        ]

    def run(self):
        counts_hdf5 = Path(self.input().path)
        run_id = Path(counts_hdf5.stem).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        r = self.cf['R']
        denoised_cr_tsv = Path(self.output()[0].path)
        standardized_cr_tsv = Path(self.output()[1].path)
        seq_dict = Path(self.seq_dict_path)
        run_dir = denoised_cr_tsv.parent
        self.setup_shell(
            run_id=run_id, commands=[gatk, r], cwd=run_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} DenoiseReadCounts'
                + f' --input {counts_hdf5}'
                + f' --standardized-copy-ratios {standardized_cr_tsv}'
                + f' --denoised-copy-ratios {denoised_cr_tsv}'
            ),
            input_files_or_dirs=counts_hdf5,
            output_files_or_dirs=[standardized_cr_tsv, denoised_cr_tsv]
        )
        if self.create_plots:
            self.run_shell(
                args=(
                    f'set -e && {gatk} PlotDenoisedCopyRatios'
                    + f' --standardized-copy-ratios {standardized_cr_tsv}'
                    + f' --denoised-copy-ratios {denoised_cr_tsv}'
                    + f' --sequence-dictionary {seq_dict}'
                    + f' --output {run_dir}'
                    + f' --output-prefix {run_id}'

                ),
                input_files_or_dirs=[
                    standardized_cr_tsv, denoised_cr_tsv, seq_dict
                ],
                output_files_or_dirs=run_dir.joinpath(f'{run_id}.denoised.png')
            )


@requires(DenoiseReadCounts)
class ModelSegments(VclineTask):
    normal_allelic_counts_tsv_path = luigi.Parameter()
    seq_dict_path = luigi.Parameter()
    case_allelic_counts_tsv_path = luigi.Parameter(default='')
    cf = luigi.DictParameter()
    dest_dir_path = luigi.Parameter(default='.')
    min_total_allele_counts = luigi.DictParameter(
        default={'case': 0, 'normal': 30}
    )
    create_plots = luigi.BoolParameter(default=True)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        output_path_prefix = str(
            Path(self.dest_dir_path).joinpath(
                Path(
                    self.create_matched_id(
                        self.case_allelic_counts_tsv_path,
                        self.normal_allelic_counts_tsv_path
                    ) if self.case_allelic_counts_tsv_path else
                    Path(self.normal_allelic_counts_tsv_path).stem
                ).stem
            )
        )
        return [
            luigi.LocalTarget(f'{output_path_prefix}.{s}')
            for s in ['cr.seg', 'hets.tsv', 'modelFinal.seg']
        ]

    def run(self):
        output_files = [Path(o.path) for o in self.output()]
        run_id = Path(output_files[0].stem).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        r = self.cf['R']
        denoised_cr_tsv = Path(self.input()[0].path)
        normal_allelic_counts_tsv = Path(
            self.normal_allelic_counts_tsv_path
        ).resolve()
        run_dir = output_files[0].parent
        if self.case_allelic_counts_tsv_path:
            case_allelic_counts_tsv = Path(
                self.case_allelic_counts_tsv_path
            ).resolve()
            allelic_count_args = (
                f' --allelic-counts {case_allelic_counts_tsv}'
                + f' --normal-allelic-counts {normal_allelic_counts_tsv}'
                + ''.join(
                    f' --minimum-total-allele-count-{k} {v}'
                    for k, v in self.min_total_allele_counts.items()
                )
            )
            input_files = [
                denoised_cr_tsv, case_allelic_counts_tsv,
                normal_allelic_counts_tsv
            ]
        else:
            allelic_count_args = (
                f' --allelic-counts {normal_allelic_counts_tsv}'
                + ' --minimum-total-allele-count-case {}'.format(
                    self.min_total_allele_counts['normal']
                )
            )
            input_files = [denoised_cr_tsv, normal_allelic_counts_tsv]
        self.setup_shell(
            run_id=run_id, commands=[gatk, r], cwd=run_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }

        )
        self.run_shell(
            args=(
                f'set -e && {gatk} ModelSegments'
                + f' --denoised-copy-ratios {denoised_cr_tsv}'
                + allelic_count_args
                + f' --output-prefix {run_id}'
                + f' --output {run_dir}'
            ),
            input_files_or_dirs=input_files,
            output_files_or_dirs=[*output_files, run_dir]
        )
        if self.create_plots:
            plots_dir = run_dir.joinpath(f'{run_id}.plots')
            seq_dict = Path(self.seq_dict_path)
            het_allelic_counts_tsv = output_files[1]
            modeled_segments = output_files[2]
            self.run_shell(
                args=(
                    f'set -e && {gatk} PlotModeledSegments'
                    + f' --denoised-copy-ratios {denoised_cr_tsv}'
                    + f' --allelic-counts {het_allelic_counts_tsv}'
                    + f' --segments {modeled_segments}'
                    + f' --sequence-dictionary {seq_dict}'
                    + f' --output {plots_dir}'
                    + f' --output-prefix {run_id}'
                ),
                input_files_or_dirs=[
                    denoised_cr_tsv, het_allelic_counts_tsv, modeled_segments,
                    seq_dict
                ],
                output_files_or_dirs=[
                    plots_dir.joinpath(f'{run_id}.modeled.png'), plots_dir
                ]
            )


@requires(ModelSegments)
class CallCopyRatioSegments(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.seg$', '.called.seg', self.input()[0].path)
        )

    def run(self):
        cr_seg = Path(self.input()[0].path)
        run_id = cr_seg.stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        gatk = self.cf['gatk']
        output_seg = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=output_seg.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} CallCopyRatioSegments'
                + f' --input {cr_seg} --output {output_seg}'
            ),
            input_files_or_dirs=cr_seg, output_files_or_dirs=output_seg
        )


@requires(PrepareCramTumor, PrepareCramNormal, PreprocessIntervals,
          FetchReferenceFasta, CreateSequenceDictionary,
          CollectAllelicCountsTumor, CollectAllelicCountsNormal)
class CallCopyRatioSegmentsTumor(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
            self.create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return luigi.LocalTarget(
            run_dir.joinpath(f'{run_dir.name}.cr.called.seg')
        )

    def run(self):
        yield CallCopyRatioSegments(
            cram_path=self.input()[0][0].path,
            preprocessed_interval_path=self.input()[2].path,
            fa_path=self.input()[3][0].path,
            seq_dict_path=self.input()[4].path,
            case_allelic_counts_tsv_path=self.input()[5].path,
            normal_allelic_counts_tsv_path=self.input()[6].path,
            dest_dir_path=str(Path(self.output().path).parent), cf=self.cf,
            n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            sh_config=self.sh_config
        )


@requires(PrepareCramTumor, PrepareCramNormal, PreprocessIntervals,
          FetchReferenceFasta, CreateSequenceDictionary,
          CollectAllelicCountsNormal)
class CallCopyRatioSegmentsNormal(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        input_cram_paths = [i[0].path for i in self.input()[0:2]]
        return luigi.LocalTarget(
            Path(self.cf['somatic_cnv_gatk_dir_path']).joinpath(
                self.create_matched_id(*input_cram_paths)
            ).joinpath(
                Path(input_cram_paths[1]).stem + '.cr.called.seg'
            )
        )

    def run(self):
        yield CallCopyRatioSegments(
            cram_path=self.input()[1][0].path,
            preprocessed_interval_path=self.input()[2].path,
            fa_path=self.input()[3][0].path,
            seq_dict_path=self.input()[4].path,
            normal_allelic_counts_tsv_path=self.input()[5].path,
            dest_dir_path=str(Path(self.output().path).parent), cf=self.cf,
            n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            sh_config=self.sh_config
        )


@requires(CallCopyRatioSegmentsTumor, CallCopyRatioSegmentsNormal)
class CallCopyRatioSegmentsMatched(luigi.WrapperTask):
    priority = 30

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
