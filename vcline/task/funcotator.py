#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask
from .haplotypecaller import FilterVariantTranches
from .mutect2 import FilterMutectCalls
from .ref import ExtractTarFile, FetchReferenceFASTA


class AnnotateGatkVCF(ShellTask):
    variant_caller = luigi.Parameter(default='mutect2')
    funcotator_data_source_tar_path = luigi.Parameter()
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    known_indel_vcf_paths = luigi.ListParameter()
    hapmap_vcf_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def requires(self):
        assert self.variant_caller in {'mutect2', 'haplotypecaller'}
        return [
            (
                FilterMutectCalls(
                    fq_list=self.fq_list, read_groups=self.read_groups,
                    sample_names=self.sample_names,
                    ref_fa_paths=self.ref_fa_paths,
                    dbsnp_vcf_path=self.dbsnp_vcf_path,
                    known_indel_vcf_paths=self.known_indel_vcf_paths,
                    gnomad_vcf_path=self.gnomad_vcf_path, cf=self.cf
                ) if self.variant_caller == 'mutect2'
                else FilterVariantTranches(
                    fq_list=self.fq_list, read_groups=self.read_groups,
                    sample_names=self.sample_names,
                    ref_fa_paths=self.ref_fa_paths,
                    dbsnp_vcf_path=self.dbsnp_vcf_path,
                    known_indel_vcf_paths=self.known_indel_vcf_paths,
                    hapmap_vcf_path=self.hapmap_vcf_path, cf=self.cf
                )
            ),
            FetchReferenceFASTA(ref_fa_paths=self.ref_fa_paths, cf=self.cf),
            ExtractTarFile(
                tar_path=self.funcotator_data_source_tar_path,
                ref_dir_path=self.cf['ref_dir_path'],
                log_dir_path=self.cf['log_dir_path']
            )
        ]

    def output(self):
        return [
            luigi.LocalTarget(
                re.sub(
                    r'(\.vcf|\.vcf\.gz)$', f'.funcotator.{s}',
                    self.input()[0][0].path
                )
            ) for s in ['vcf.gz', 'vcf.gz.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Annotate variants with Funcotator:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        input_vcf_path = self.input()[0][0].path
        fa_path = self.input()[1].path
        data_src_dir_path = self.input()[2].path
        ref_version = self.cf['ref_version']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=str(Path(output_vcf_path).parent)
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} Funcotator'
                + f' --variant {input_vcf_path}'
                + f' --reference {fa_path}'
                + f' --ref-version {ref_version}'
                + f' --data-sources-path {data_src_dir_path}'
                + f' --output {output_vcf_path}'
                + f' --output-file-format VCF'
            ),
            input_files=[input_vcf_path, fa_path, data_src_dir_path],
            output_files=[output_vcf_path, f'{output_vcf_path}.tbi']
        )


if __name__ == '__main__':
    luigi.run()
