#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from ftarc.task.base import ShellTask
from ftarc.task.picard import CreateSequenceDictionary
from ftarc.task.resource import FetchReferenceFASTA

from .bcftools import NormalizeVCF


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
                tar_path=self.data_src_tar_path,
                dest_dir_path=str(
                    Path(self.cf['postproc_dir_path']).joinpath('data_sources')
                ),
                cf=self.cf
            ),
            FetchReferenceFASTA(ref_fa_path=self.ref_fa_path, cf=self.cf),
            CreateSequenceDictionary(ref_fa_path=self.ref_fa_path, cf=self.cf)
        ]

    def output(self):
        is_seg = (self.output_file_format == 'SEG')
        output_vcf = Path(self.cf['postproc_dir_path']).joinpath(
            'annotation'
        ).joinpath(
            re.sub(r'\.vcf$', '', Path(self.input_vcf_path).stem)
            + ('.norm' if self.normalize_vcf else '') + '.funcotator.'
            + self.output_file_format.lower()
            + ('.tsv' if is_seg else '.gz')
        )
        return [
            luigi.LocalTarget(f'{output_vcf}{s}')
            for s in ['', ('.gene_list.txt' if is_seg else '.tbi')]
        ]

    def run(self):
        dest_dir = Path(self.output()[0].path).parent
        yield Funcotator(
            input_vcf_path=self.input_vcf_path,
            data_src_dir_path=self.input()[0].path,
            fa_path=self.input()[1][0].path,
            ref_version=self.cf['ref_version'], dest_dir_path=str(dest_dir),
            normalize_vcf=self.normalize_vcf,
            norm_dir_path=str(dest_dir.parent.joinpath('normalization')),
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


class ExtractTarFile(ShellTask):
    tar_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter()
    cf = luigi.DictParameter()
    recursive = luigi.BoolParameter(default=True)
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).joinpath(
                Path(Path(self.tar_path).stem).stem
            )
        )

    def run(self):
        target_dir = Path(self.output().path)
        run_id = target_dir.name
        self.print_log(f'Extract a tar file:\t{run_id}')
        pigz = self.cf['pigz']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[pigz, pbzip2], cwd=self.dest_dir_path
        )
        self._tar_xf(
            tar_path=self.tar_path, pigz=pigz, pbzip2=pbzip2, n_cpu=n_cpu,
            input_files_or_dirs=self.tar_path, output_files_or_dirs=target_dir
        )
        if self.recursive and target_dir.is_dir():
            for f in target_dir.iterdir():
                if f.name.endswith(('.tar.gz', '.tar.bz2')):
                    self._tar_xf(
                        tar_path=f, pigz=pigz, pbzip2=pbzip2, n_cpu=n_cpu,
                        cwd=target_dir, input_files_or_dirs=f,
                        output_files_or_dirs=target_dir.joinpath(
                            Path(Path(f).stem).stem
                        )
                    )

    def _tar_xf(self, tar_path, pigz, pbzip2, n_cpu, **kwargs):
        self.run_shell(
            args=(
                'set -eo pipefail && ' + (
                    f'{pbzip2} -p{n_cpu} -dc {tar_path}'
                    if str(tar_path).endswith('.bz2') else
                    f'{pigz} -p {n_cpu} -dc {tar_path}'
                ) + ' | tar xmvf -'
            ),
            **kwargs
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
        output_vcf = Path(self.dest_dir_path).joinpath(
            re.sub(
                r'\.vcf$', '',
                Path(
                    self.input()[0].path if self.normalize_vcf
                    else self.input_vcf_path
                ).stem
            ) + '.funcotator.' + self.output_file_format.lower()
            + ('.tsv' if self.output_file_format == 'SEG' else '.gz')
        )
        return [
            luigi.LocalTarget(f'{output_vcf}{s}') for s in [
                '',
                (
                    '.gene_list.txt' if self.output_file_format == 'SEG'
                    else '.tbi'
                )
            ]
        ]

    def run(self):
        input_vcf = Path(
            self.input()[0].path if self.normalize_vcf else self.input_vcf_path
        )
        run_id = Path(input_vcf.stem).stem
        self.print_log(f'Annotate variants with Funcotator:\t{run_id}')
        gatk_opts = (
            f' --java-options "{self.gatk_java_options}"'
            if self.gatk_java_options else ''
        )
        output_file_paths = [o.path for o in self.output()]
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=self.gatk, cwd=self.dest_dir_path,
            remove_if_failed=self.remove_if_failed, quiet=self.quiet
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk}{gatk_opts} Funcotator'
                + f' --variant {input_vcf}'
                + f' --data-sources-path {self.data_src_dir_path}'
                + f' --reference {self.fa_path}'
                + f' --ref-version {self.ref_version}'
                + f' --output {output_file_paths[0]}'
                + f' --output-file-format {self.output_file_format}'
                + (
                    ' --disable-sequence-dictionary-validation'
                    if self.disable_vcf_validation else ''
                )
            ),
            input_files_or_dirs=[
                input_vcf, self.fa_path, self.data_src_dir_path
            ],
            output_files_or_dirs=output_file_paths,
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
                tar_path=self.data_src_tar_path,
                dest_dir_path=str(
                    Path(self.cf['postproc_dir_path']).joinpath('data_sources')
                ),
                cf=self.cf
            ),
            FetchReferenceFASTA(ref_fa_path=self.ref_fa_path, cf=self.cf),
            CreateSequenceDictionary(
                ref_fa_path=self.ref_fa_path, cf=self.cf
            )
        ]

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['postproc_dir_path']).joinpath('annotation').joinpath(
                Path(self.input_seg_path).name + '.funcotator.tsv'
            )
        )

    def run(self):
        output_tsv = Path(self.output().path)
        run_id = output_tsv.stem
        self.print_log(f'Annotate segments with FuncotateSegments:\t{run_id}')
        data_src_dir_path = self.input()[0].path
        fa_path = self.input()[1][0].path
        ref_version = self.cf['ref_version']
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=str(output_tsv.parent),
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
                + f' --output {output_tsv}'
                + ' --output-file-format SEG'
            ),
            input_files_or_dirs=[
                self.input_seg_path, fa_path, data_src_dir_path
            ],
            output_files_or_dirs=output_tsv
        )


if __name__ == '__main__':
    luigi.run()
