#!/usr/bin/env python

import re
from itertools import product
from pathlib import Path

import luigi

from .base import ShellTask
from .haplotypecaller import FilterVariantTranches
from .manta import CallStructualVariantsWithManta
from .mutect2 import FilterMutectCalls
from .ref import (ExtractTarFile, FetchEvaluationIntervalList,
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
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def requires(self):
        if self.variant_caller == 'haplotypecaller':
            variant_calling = FilterVariantTranches(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                hapmap_vcf_path=self.hapmap_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.variant_caller == 'mutect2':
            variant_calling = FilterMutectCalls(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                gnomad_vcf_path=self.gnomad_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.variant_caller in {'manta_somatic', 'manta_diploid'}:
            variant_calling = CallStructualVariantsWithManta(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.variant_caller == 'strelka':
            variant_calling = CallVariantsWithStrelka(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
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
            FetchEvaluationIntervalList(
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        ]

    def output(self):
        return [
            luigi.LocalTarget(v + i) for v, i in product(
                [
                    re.sub(r'\.vcf\.gz$', '.Funcotator.vcf.gz', p)
                    for p in self._generate_input_vcf_paths()
                ],
                ['', '.tbi']
            )
        ]

    def _generate_input_vcf_paths(self):
        prefix = (
            self.variant_caller.split('_')[1]
            if '_' in self.variant_caller else None
        )
        for i in self.input()[0]:
            p = i.path
            if ((prefix is None or Path(p).name.startswith(prefix))
                    and p.endswith('.vcf.gz')):
                yield p

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
        asynchronous = (self.cf['n_cpu_per_worker'] > 1)
        input_vcf_paths = list(self._generate_input_vcf_paths())
        fa_path = self.input()[1].path
        data_src_dir_path = self.input()[2].path
        evaluation_interval_path = self.input()[3].path
        ref_version = self.cf['ref_version']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=str(Path(output_vcf_paths[0]).parent)
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {gatk}{gatk_opts} Funcotator'
                    + f' --variant {i}'
                    + f' --reference {fa_path}'
                    + f' --intervals {evaluation_interval_path}'
                    + f' --ref-version {ref_version}'
                    + f' --data-sources-path {data_src_dir_path}'
                    + f' --output {o}'
                    + ' --output-file-format VCF'
                ) for i, o in zip(input_vcf_paths, output_vcf_paths)
            ],
            input_files=[
                *input_vcf_paths, fa_path, data_src_dir_path,
                evaluation_interval_path
            ],
            output_files=[
                *output_vcf_paths, *[f'{o}.tbi' for o in output_vcf_paths]
            ],
            asynchronous=asynchronous
        )


if __name__ == '__main__':
    luigi.run()
