#!/usr/bin/env python

import re
from itertools import product
from pathlib import Path

import luigi

from .base import ShellTask
from .haplotypecaller import FilterVariantTranches
from .manta import CallStructualVariantsWithManta
from .mutect2 import FilterMutectCalls
from .ref import (CreateEvaluationIntervalList, ExtractTarFile,
                  FetchReferenceFASTA)
from .strelka import CallVariantsWithStrelka


class AnnotateVariantsWithFuncotator(ShellTask):
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
        if self.variant_caller == 'haplotypecaller':
            variant_calling = FilterVariantTranches(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                hapmap_vcf_path=self.hapmap_vcf_path, cf=self.cf
            )
        elif self.variant_caller == 'mutect2':
            variant_calling = FilterMutectCalls(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                gnomad_vcf_path=self.gnomad_vcf_path, cf=self.cf
            )
        elif self.variant_caller in {'manta_somatic', 'manta_diploid'}:
            variant_calling = CallStructualVariantsWithManta(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths, cf=self.cf
            )
        elif self.variant_caller == 'strelka':
            variant_calling = CallVariantsWithStrelka(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths, cf=self.cf
            )
        else:
            raise ValueError(f'invalid variant_caller: {self.variant_caller}')
        return [
            variant_calling,
            FetchReferenceFASTA(ref_fa_paths=self.ref_fa_paths, cf=self.cf),
            ExtractTarFile(
                tar_path=self.funcotator_data_source_tar_path,
                ref_dir_path=self.cf['ref_dir_path'],
                log_dir_path=self.cf['log_dir_path']
            ),
            CreateEvaluationIntervalList(
                ref_fa_paths=self.ref_fa_paths, cf=self.cf
            )
        ]

    def output(self):
        if self.variant_caller.startswith('manta_'):
            input_vcf_paths = [
                i.path for i in self.input()[0] if
                Path(i.path).name.startswith(self.variant_caller.split('_')[1])
            ]
        else:
            input_vcf_paths = [
                i.path for i in self.input()[0] if i.path.endswith('.vcf.gz')
            ]
        return [
            luigi.LocalTarget(v + i) for v, i in product(
                [
                    re.sub(r'\.vcf\.gz$', '.Funcotator.vcf.gz', p)
                    for p in input_vcf_paths
                ],
                ['', '.tbi']
            )
        ]

    def run(self):
        output_vcf_paths = [
            o.path for o in self.output() if o.path.endswith('.vcf.gz')
        ]
        run_id = (
            Path(output_vcf_paths[0]).parent.parent.parent.name
            if self.variant_caller in {'manta', 'strelka'} else
            '.'.join(Path(output_vcf_paths[0]).name.split('.')[:-3])
        )
        self.print_log(f'Annotate variants with Funcotator:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        fa_path = self.input()[1].path
        data_src_dir_path = self.input()[2].path
        evaluation_interval_path = self.input()[3].path
        ref_version = self.cf['ref_version']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=str(Path(output_vcf_paths[0]).parent)
        )
        for output_vcf_path in output_vcf_paths:
            input_vcf_path = re.sub(
                r'\.Funcotator\.vcf\.gz', '.vcf.gz', output_vcf_path
            )
            self.run_shell(
                args=(
                    f'set -e && {gatk}{gatk_opts} Funcotator'
                    + f' --variant {input_vcf_path}'
                    + f' --reference {fa_path}'
                    + f' --intervals {evaluation_interval_path}'
                    + f' --ref-version {ref_version}'
                    + f' --data-sources-path {data_src_dir_path}'
                    + f' --output {output_vcf_path}'
                    + f' --output-file-format VCF'
                ),
                input_files=[
                    input_vcf_path, fa_path, data_src_dir_path,
                    evaluation_interval_path
                ],
                output_files=[output_vcf_path, f'{output_vcf_path}.tbi']
            )


if __name__ == '__main__':
    luigi.run()
