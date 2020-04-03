#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask
from .bcftools import NormalizeVCF
from .ref import (CreateFASTAIndex, CreateSequenceDictionary, ExtractTarFile,
                  FetchReferenceFASTA)


class AnnotateVariantsWithFuncotator(luigi.WrapperTask):
    input_vcf_path = luigi.Parameter()
    data_src_tar_path = luigi.Parameter()
    ref_fa_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = 10

    def requires(self):
        return [
            ExtractTarFile(
                tar_path=self.data_src_tar_path, cf=self.cf
            ),
            FetchReferenceFASTA(ref_fa_paths=self.ref_fa_paths, cf=self.cf),
            CreateFASTAIndex(ref_fa_paths=self.ref_fa_paths, cf=self.cf),
            CreateSequenceDictionary(
                ref_fa_paths=self.ref_fa_paths, cf=self.cf
            )
        ]

    def output(self):
        return RunFuncotator(
            input_vcf_path=self.input_vcf_path,
            data_src_dir_path=self.input()[0].path,
            fa_path=self.input()[1].path, ref_version=self.cf['ref_version'],
            dest_dir_path=self.cf['funcotator_dir_path'],
            normalize_vcf=self.normalize_vcf,
            norm_dir_path=self.cf['bcftools_dir_path'],
            n_cpu=self.cf['n_cpu_per_worker'],
            memory_mb=self.cf['memory_mb_per_worker'], gatk=self.cf['gatk'],
            gatk_java_options=self.cf['gatk_java_options'],
            bcftools=self.cf['bcftools'], log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            disable_vcf_validation=False
        ).output()


class RunFuncotator(ShellTask):
    input_vcf_path = luigi.Parameter()
    data_src_dir_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    ref_version = luigi.Parameter(default='hg38')
    dest_dir_path = luigi.Parameter(default='.')
    normalize_vcf = luigi.BoolParameter(default=False)
    norm_dir_path = luigi.Parameter(default='')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.IntParameter(default=(4 * 1024))
    gatk = luigi.Parameter()
    gatk_java_options = luigi.Parameter(default='')
    bcftools = luigi.Parameter()
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    disable_vcf_validation = luigi.BoolParameter(default=False)
    priority = 10

    def requires(self):
        if self.normalize_vcf:
            return NormalizeVCF(
                input_vcf_path=self.input_vcf_path, fa_path=self.fa_path,
                dest_dir_path=(self.norm_dir_path or self.dest_dir_path),
                n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                bcftools=self.bcftools,
                log_dir_path=(self.log_dir_path or None),
                remove_if_failed=self.remove_if_failed
            )
        else:
            return super().requires()

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.dest_dir_path).joinpath(
                        re.sub(
                            r'\.vcf$', '',
                            Path(
                                self.input()[0].path if self.normalize_vcf
                                else self.input_vcf_path
                            ).stem
                        ) + f'.funcotator.vcf.gz{s}'
                    )
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Annotate variants with Funcotator:\t{run_id}')
        input_vcf_path = (
            self.input()[0].path if self.normalize_vcf else self.input_vcf_path
        )
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.gatk, cwd=self.dest_dir_path,
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=(
                f'set -e && '
                + (
                    f'{self.gatk} --java-options "{self.gatk_java_options}"'
                    if self.gatk_java_options else self.gatk
                ) + ' Funcotator'
                + f' --variant {input_vcf_path}'
                + f' --reference {self.fa_path}'
                + f' --ref-version {self.ref_version}'
                + f' --data-sources-path {self.data_src_dir_path}'
                + f' --output {output_vcf_path}'
                + ' --output-file-format VCF'
                + (
                    ' --disable-sequence-dictionary-validation'
                    if self.disable_vcf_validation else ''
                )
            ),
            input_files_or_dirs=[
                input_vcf_path, self.fa_path, self.data_src_dir_path
            ],
            output_files_or_dirs=[o.path for o in self.output()],
        )


if __name__ == '__main__':
    luigi.run()
