#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import ShellTask
from .bcftools import NormalizeVCF
from .ref import CreateSequenceDictionary, FetchReferenceFASTA
from .samtools import tabix_tbi


class AnnotateVariantsWithSnpEff(luigi.Task):
    input_vcf_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    snpeff_config_path = luigi.Parameter()
    cf = luigi.DictParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = 10

    def requires(self):
        return [
            FetchReferenceFASTA(ref_fa_path=self.ref_fa_path, cf=self.cf),
            CreateSequenceDictionary(ref_fa_path=self.ref_fa_path, cf=self.cf)
        ]

    def output(self):
        output_vcf = Path(self.cf['postproc_dir_path']).joinpath(
            'annotation'
        ).joinpath(
            re.sub(r'\.vcf$', '', Path(self.input_vcf_path).stem)
            + ('.norm' if self.normalize_vcf else '') + '.snpeff.vcf.gz'
        )
        return [luigi.LocalTarget(f'{output_vcf}{s}') for s in ['', '.tbi']]

    def run(self):
        dest_dir = Path(self.output()[0].path).parent
        yield SnpEff(
            input_vcf_path=self.input_vcf_path,
            fa_path=self.input()[0][0].path,
            snpeff_config_path=self.snpeff_config_path,
            ref_version=self.cf['ref_version'], dest_dir_path=str(dest_dir),
            normalize_vcf=self.normalize_vcf,
            norm_dir_path=str(dest_dir.parent.joinpath('normalization')),
            bcftools=self.cf['bcftools'], snpeff=self.cf['snpEff'],
            bgzip=self.cf['bgzip'], tabix=self.cf['tabix'],
            n_cpu=self.cf['n_cpu_per_worker'],
            memory_mb=self.cf['memory_mb_per_worker'],
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )


class SnpEff(ShellTask):
    input_vcf_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    snpeff_config_path = luigi.Parameter()
    ref_version = luigi.Parameter(default='hg38')
    snpeff_genome_version = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    normalize_vcf = luigi.BoolParameter(default=False)
    norm_dir_path = luigi.Parameter(default='')
    bcftools = luigi.Parameter()
    snpeff = luigi.Parameter()
    bgzip = luigi.Parameter()
    tabix = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.IntParameter(default=(4 * 1024))
    log_dir_path = luigi.Parameter(default='')
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
            ) + '.snpeff.vcf.gz'
        )
        return [luigi.LocalTarget(f'{output_vcf}{s}') for s in ['', '.tbi']]

    def run(self):
        input_vcf = Path(
            self.input()[0].path if self.normalize_vcf else self.input_vcf_path
        )
        run_id = Path(input_vcf.stem).stem
        self.print_log(f'Annotate variants with SnpEff:\t{run_id}')
        snpeff_config = Path(self.snpeff_config_path).resolve()
        tmp_dir = Path(self.dest_dir_path).joinpath(run_id)
        tmp_file_paths = [
            tmp_dir.joinpath(n) for n
            in ['snpeff.vcf.gz', 'snpEff_genes.txt', 'snpEff_summary.html']
        ]
        genome_version = (
            self.snpeff_genome_version or [
                o.name
                for o in snpeff_config.parent.joinpath('data').iterdir()
                if (
                    o.name.startswith(
                        {'hg38': 'GRCh38', 'hg19': 'GRCh38'}[self.ref_version]
                    ) and o.is_dir()
                )
            ][0]
        )
        output_vcf_path = self.output()[0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=[self.snpeff, self.bgzip, self.tabix],
            cwd=self.dest_dir_path, remove_if_failed=self.remove_if_failed,
            quiet=self.quiet
        )
        self.run_shell(args=f'mkdir {tmp_dir}', output_files_or_dirs=tmp_dir)
        self.run_shell(
            args=(
                f'set -e && cd {tmp_dir} && '
                + f'{self.snpeff} -verbose -config {snpeff_config}'
                + f' {genome_version} {input_vcf}'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {tmp_file_paths[0]}'
            ),
            input_files_or_dirs=[input_vcf, tmp_dir],
            output_files_or_dirs=[tmp_file_paths[0], tmp_dir]
        )
        for p in tmp_file_paths:
            i = Path(p)
            if i.is_file():
                o = Path(self.dest_dir_path).joinpath(f'{run_id}.{i.name}')
                self.run_shell(
                    args=f'mv {i} {o}', input_files_or_dirs=str(i),
                    output_files_or_dirs=str(o)
                )
        self.run_shell(args=f'rm -rf {tmp_dir}', input_files_or_dirs=tmp_dir)
        tabix_tbi(
            shelltask=self, tabix=self.tabix, tsv_path=output_vcf_path,
            preset='vcf'
        )


if __name__ == '__main__':
    luigi.run()
