#!/usr/bin/env python

from itertools import product
from pathlib import Path

import luigi
from luigi.util import requires

from .base import ShellTask
from .delly import CallStructualVariantsWithDelly
from .haplotypecaller import FilterVariantTranches
from .lumpy import CallStructualVariantsWithLumpy
from .manta import CallStructualVariantsWithManta
from .mutect2 import FilterMutectCalls
from .ref import (ExtractFuncotatorTarFile, FetchEvaluationIntervalList,
                  FetchReferenceFASTA)
from .strelka import (CallGermlineVariantsWithStrelka,
                      CallSomaticVariantsWithStrelka)


class RunVariantCaller(luigi.WrapperTask):
    variant_caller = luigi.Parameter()
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    hapmap_vcf_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def requires(self):
        if self.variant_caller == 'haplotypecaller':
            return FilterVariantTranches(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                hapmap_vcf_path=self.hapmap_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.variant_caller == 'mutect2':
            return FilterMutectCalls(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                gnomad_vcf_path=self.gnomad_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.variant_caller == 'manta_somatic':
            return CallStructualVariantsWithManta(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.variant_caller == 'strelka_somatic':
            return CallSomaticVariantsWithStrelka(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.variant_caller == 'strelka_variants':
            return CallGermlineVariantsWithStrelka(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif self.variant_caller == 'delly':
            return CallStructualVariantsWithDelly(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                cf=self.cf
            )
        elif self.variant_caller == 'lumpy':
            return CallStructualVariantsWithLumpy(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                cf=self.cf
            )
        else:
            raise ValueError(f'invalid variant_caller: {self.variant_caller}')

    def output(self):
        return self.input()


@requires(RunVariantCaller, FetchReferenceFASTA, ExtractFuncotatorTarFile,
          FetchEvaluationIntervalList)
class AnnotateVariantsWithFuncotator(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(v + i) for v, i
            in product(self._generate_output_vcf_paths(), ['', '.tbi'])
        ]

    def _generate_input_vcf_paths(self):
        prefix = (
            self.variant_caller.split('_')[1]
            if '_' in self.variant_caller else None
        )
        for i in self.input()[0]:
            p = i.path
            if (p.endswith('.vcf.gz')
                    and (prefix is None or Path(p).name.startswith(prefix))):
                yield p

    def _generate_output_vcf_paths(self):
        vc = self.variant_caller.split('_')[0]
        for p in self._generate_input_vcf_paths():
            yield str(
                Path(self.cf['funcotator_dir_path']).joinpath(
                    '.'.join(
                        (
                            [Path(p).parent.parent.parent.name, vc]
                            if vc in {'manta', 'strelka'} else list()
                        ) + [Path(Path(p).stem).stem, 'funcotator.vcf.gz']
                    )
                )
            )

    def run(self):
        output_vcf_paths = list(self._generate_output_vcf_paths())
        run_id = '.'.join(Path(output_vcf_paths[0]).name.split('.')[:-3])
        self.print_log(f'Annotate variants with Funcotator:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        input_vcf_paths = list(self._generate_input_vcf_paths())
        fa_path = self.input()[1].path
        data_src_dir_path = self.input()[2].path
        evaluation_interval_path = self.input()[3].path
        ref_version = self.cf['ref_version']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['funcotator_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
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
            input_files_or_dirs=[
                *input_vcf_paths, fa_path, data_src_dir_path,
                evaluation_interval_path
            ],
            output_files_or_dirs=[
                *output_vcf_paths, *[f'{o}.tbi' for o in output_vcf_paths]
            ],
            asynchronous=(
                1 < len(input_vcf_paths) <= self.cf['n_cpu_per_worker']
            )
        )


if __name__ == '__main__':
    luigi.run()
