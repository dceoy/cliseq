#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask
from .bcftools import NormalizeVCF
from .ref import CreateSequenceDictionary, ExtractTarFile, FetchReferenceFASTA


class FuncotateVariants(luigi.Task):
    input_vcf_path = luigi.Parameter()
    data_src_tar_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    output_file_format = luigi.Parameter(default='VCF')
    priority = 10

    def requires(self):
        return [
            ExtractTarFile(
                tar_path=self.data_src_tar_path, cf=self.cf
            ),
            FetchReferenceFASTA(ref_fa_path=self.ref_fa_path, cf=self.cf),
            CreateSequenceDictionary(
                ref_fa_path=self.ref_fa_path, cf=self.cf
            )
        ]

    def output(self):
        suffixes = [
            '{0}.funcotator.{1}.{2}'.format(
                ('.norm' if self.normalize_vcf else ''),
                self.output_file_format.lower(),
                (f'tsv{s}' if self.output_file_format == 'SEG' else f'gz{s}')
            ) for s in [
                '',
                (
                    '.gene_list.txt' if self.output_file_format == 'SEG'
                    else '.tbi'
                )
            ]
        ]
        return [
            luigi.LocalTarget(
                Path(self.cf['postproc_funcotator_dir_path']).joinpath(
                    re.sub(r'\.vcf$', s, Path(self.input_vcf_path).stem)
                )
            ) for s in suffixes
        ]

    def run(self):
        yield Funcotator(
            input_vcf_path=self.input_vcf_path,
            data_src_dir_path=self.input()[0].path,
            fa_path=self.input()[1][0].path,
            ref_version=self.cf['ref_version'],
            dest_dir_path=self.cf['postproc_funcotator_dir_path'],
            normalize_vcf=self.normalize_vcf,
            norm_dir_path=self.cf['postproc_bcftools_dir_path'],
            bcftools=self.cf['bcftools'], gatk=self.cf['gatk'],
            gatk_java_options=self.cf['gatk_java_options'],
            n_cpu=self.cf['n_cpu_per_worker'],
            memory_mb=self.cf['memory_mb_per_worker'],
            log_dir_path=self.cf['log_dir_path'],
            output_file_format=self.output_file_format,
            disable_vcf_validation=False,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )


class Funcotator(ShellTask):
    input_vcf_path = luigi.Parameter()
    data_src_dir_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    ref_version = luigi.Parameter(default='hg38')
    dest_dir_path = luigi.Parameter(default='.')
    normalize_vcf = luigi.BoolParameter(default=False)
    norm_dir_path = luigi.Parameter(default='')
    bcftools = luigi.Parameter()
    gatk = luigi.Parameter()
    gatk_java_options = luigi.Parameter(default='')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.IntParameter(default=(4 * 1024))
    log_dir_path = luigi.Parameter(default='')
    output_file_format = luigi.Parameter(default='VCF')
    disable_vcf_validation = luigi.BoolParameter(default=False)
    remove_if_failed = luigi.BoolParameter(default=True)
    quiet = luigi.BoolParameter(default=False)
    priority = 10

    def requires(self):
        if self.normalize_vcf:
            return NormalizeVCF(
                input_vcf_path=self.input_vcf_path, fa_path=self.fa_path,
                dest_dir_path=(self.norm_dir_path or self.dest_dir_path),
                n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                bcftools=self.bcftools,
                log_dir_path=(self.log_dir_path or None),
                remove_if_failed=self.remove_if_failed, quiet=self.quiet
            )
        else:
            return super().requires()

    def output(self):
        suffixes = [
            '.funcotator.{0}.{1}'.format(
                self.output_file_format.lower(),
                (f'tsv{s}' if self.output_file_format == 'SEG' else f'gz{s}')
            ) for s in [
                '',
                (
                    '.gene_list.txt' if self.output_file_format == 'SEG'
                    else '.tbi'
                )
            ]
        ]
        return [
            luigi.LocalTarget(
                Path(self.dest_dir_path).joinpath(
                    re.sub(
                        r'\.vcf$', s,
                        Path(
                            self.input()[0].path if self.normalize_vcf
                            else self.input_vcf_path
                        ).stem
                    )
                )
            ) for s in suffixes
        ]

    def run(self):
        output_file_path = self.output()[0].path
        run_id = '.'.join(Path(output_file_path).name.split('.')[:-3])
        self.print_log(f'Annotate variants with Funcotator:\t{run_id}')
        gatk_opts = (
            f' --java-options "{self.gatk_java_options}"'
            if self.gatk_java_options else ''
        )
        input_vcf_path = (
            self.input()[0].path if self.normalize_vcf else self.input_vcf_path
        )
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.gatk, cwd=self.dest_dir_path,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk}{gatk_opts} Funcotator'
                + f' --variant {input_vcf_path}'
                + f' --data-sources-path {self.data_src_dir_path}'
                + f' --reference {self.fa_path}'
                + f' --ref-version {self.ref_version}'
                + f' --output {output_file_path}'
                + f' --output-file-format {self.output_file_format}'
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


class FuncotateSegments(ShellTask):
    input_seg_path = luigi.Parameter()
    data_src_tar_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def requires(self):
        return [
            ExtractTarFile(
                tar_path=self.data_src_tar_path, cf=self.cf
            ),
            FetchReferenceFASTA(ref_fa_path=self.ref_fa_path, cf=self.cf),
            CreateSequenceDictionary(
                ref_fa_path=self.ref_fa_path, cf=self.cf
            )
        ]

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['postproc_funcotator_dir_path']).joinpath(
                Path(self.input_seg_path).name + '.funcotator.tsv'
            )
        )

    def run(self):
        run_id = Path(self.input_seg_path).stem
        self.print_log(f'Annotate segments with FuncotateSegments:\t{run_id}')
        output_tsv_path = self.output().path
        data_src_dir_path = self.input()[0].path
        fa_path = self.input()[1][0].path
        ref_version = self.cf['ref_version']
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['postproc_funcotator_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} FuncotateSegments'
                + f' --segments {self.input_seg_path}'
                + f' --data-sources-path {data_src_dir_path}'
                + f' --reference {fa_path}'
                + f' --ref-version {ref_version}'
                + f' --output {output_tsv_path}'
                + ' --output-file-format SEG'
            ),
            input_files_or_dirs=[
                self.input_seg_path, fa_path, data_src_dir_path
            ],
            output_files_or_dirs=output_tsv_path
        )


if __name__ == '__main__':
    luigi.run()
