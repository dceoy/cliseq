#!/usr/bin/env python

from pathlib import Path

import luigi

from .base import ShellTask
from .bcftools import NormalizeVCF
from .ref import (ExtractFuncotatorTarFile, FetchEvaluationIntervalList,
                  FetchReferenceFASTA)


class AnnotateVariantsWithFuncotator(ShellTask):
    input_vcf_path = luigi.Parameter()
    data_source_tar_path = luigi.Parameter()
    ref_fa_paths = luigi.ListParameter()
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = 10

    def requires(self):
        return [
            FetchReferenceFASTA(ref_fa_paths=self.ref_fa_paths, cf=self.cf),
            ExtractFuncotatorTarFile(
                tar_path=self.data_source_tar_path, cf=self.cf
            ),
            FetchEvaluationIntervalList(
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            ),
            *(
                [
                    NormalizeVCF(
                        input_vcf_path=self.input_vcf_path,
                        output_vcf_path=str(
                            Path(self.cf['funcotator_dir_path']).joinpath(
                                Path(Path(self.input_vcf_path).stem).stem
                                + ('.norm' if self.normalize_vcf else '')
                                + '.vcf.gz'
                            )
                        ),
                        ref_fa_paths=self.ref_fa_paths, cf=self.cf
                    )
                ] if self.normalize_vcf else list()
            )
        ]

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['funcotator_dir_path']).joinpath(
                    Path(Path(self.input_vcf_path).stem).stem
                    + ('.norm' if self.normalize_vcf else '')
                    + f'.funcotator.vcf.gz{s}'
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Annotate variants with Funcotator:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        fa_path = self.input()[0].path
        data_src_dir_path = self.input()[1].path
        evaluation_interval_path = self.input()[2].path
        input_vcf_path = (
            self.input()[3].path if self.normalize_vcf else self.input_vcf_path
        )
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['funcotator_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                    f'set -e && {gatk}{gatk_opts} Funcotator'
                    + f' --variant {input_vcf_path}'
                    + f' --reference {fa_path}'
                    + f' --intervals {evaluation_interval_path}'
                    + ' --ref-version {}'.format(self.cf['ref_version'])
                    + f' --data-sources-path {data_src_dir_path}'
                    + f' --output {output_vcf_path}'
                    + ' --output-file-format VCF'
            ),
            input_files_or_dirs=[
                input_vcf_path, fa_path, data_src_dir_path,
                evaluation_interval_path
            ],
            output_files_or_dirs=[o.path for o in self.output()],
        )


if __name__ == '__main__':
    luigi.run()
