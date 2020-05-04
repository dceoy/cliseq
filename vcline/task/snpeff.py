#!/usr/bin/env python

import re
from pathlib import Path

import luigi

from .base import BaseTask, ShellTask
from .bcftools import BcftoolsIndex, NormalizeVCF
from .ref import CreateSequenceDictionary, CreateSymlinks, FetchReferenceFASTA


class AnnotateVariantsWithSnpEff(BaseTask):
    input_vcf_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    snpeff_config_path = luigi.Parameter()
    cf = luigi.DictParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = 10

    def requires(self):
        return [
            FetchReferenceFASTA(ref_fa_path=self.ref_fa_path, cf=self.cf),
            CreateSequenceDictionary(ref_fa_path=self.ref_fa_path, cf=self.cf),
            CreateSymlinks(
                src_paths=[self.snpeff_config_path],
                dest_dir_path=self.cf['ref_dir_path'], cf=self.cf
            )
        ]

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['postproc_snpeff_dir_path']).joinpath(
                    re.sub(
                        r'\.vcf$',
                        '{0}.snpeff.vcf.gz{1}'.format(
                            ('.norm' if self.normalize_vcf else ''), s
                        ),
                        Path(self.input_vcf_path).stem
                    )
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        yield SnpEff(
            input_vcf_path=self.input_vcf_path,
            fa_path=self.input()[0][0].path,
            snpeff_config_path=self.input()[2][0].path,
            ref_version=self.cf['ref_version'],
            dest_dir_path=self.cf['postproc_snpeff_dir_path'],
            normalize_vcf=self.normalize_vcf,
            norm_dir_path=self.cf['postproc_bcftools_dir_path'],
            bcftools=self.cf['bcftools'], snpeff=self.cf['snpEff'],
            bgzip=self.cf['bgzip'], n_cpu=self.cf['n_cpu_per_worker'],
            memory_mb=self.cf['memory_mb_per_worker'],
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
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
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.IntParameter(default=(4 * 1024))
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
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
                Path(self.dest_dir_path).joinpath(
                    re.sub(
                        r'\.vcf$', f'.snpeff.vcf.gz{s}',
                        Path(
                            self.input()[0].path if self.normalize_vcf
                            else self.input_vcf_path
                        ).stem
                    )
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = '.'.join(Path(output_vcf_path).name.split('.')[:-3])
        self.print_log(f'Annotate variants with SnpEff:\t{run_id}')
        input_vcf_path = (
            self.input()[0].path if self.normalize_vcf else self.input_vcf_path
        )
        config_path = str(Path(self.snpeff_config_path).resolve())
        output_prefix = Path(Path(input_vcf_path).stem).stem
        tmp_dir_path = str(Path(self.dest_dir_path).joinpath(output_prefix))
        tmp_file_paths = [
            str(Path(tmp_dir_path).joinpath(n)) for n
            in ['snpeff.vcf.gz', 'snpEff_genes.txt', 'snpEff_summary.html']
        ]
        genome_version = (
            self.snpeff_genome_version or [
                o.name
                for o in Path(config_path).parent.joinpath('data').iterdir()
                if (
                    o.name.startswith(
                        {'hg38': 'GRCh38', 'hg19': 'GRCh38'}[self.ref_version]
                    ) and o.is_dir()
                )
            ][0]
        )
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=[self.snpeff, self.bgzip],
            cwd=self.dest_dir_path, remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=f'mkdir {tmp_dir_path}', output_files_or_dirs=tmp_dir_path
        )
        self.run_shell(
            args=(
                f'set -e && cd {tmp_dir_path} && '
                + f'{self.snpeff} -verbose -config {config_path}'
                + f' {genome_version} {input_vcf_path}'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {tmp_file_paths[0]}'
            ),
            input_files_or_dirs=[input_vcf_path, tmp_dir_path],
            output_files_or_dirs=[tmp_file_paths[0], tmp_dir_path]
        )
        for p in tmp_file_paths:
            i = Path(p)
            if i.is_file():
                o = Path(self.dest_dir_path).joinpath(
                    f'{output_prefix}.{i.name}'
                )
                self.run_shell(
                    args=f'mv {i} {o}', input_files_or_dirs=str(i),
                    output_files_or_dirs=str(o)
                )
        self.run_shell(
            args=f'rm -rf {tmp_dir_path}', input_files_or_dirs=tmp_dir_path
        )
        yield BcftoolsIndex(
            vcf_path=output_vcf_path, bcftools=self.bcftools,
            n_cpu=self.n_cpu, tbi=True, log_dir_path=self.log_dir_path,
            remove_if_failed=self.remove_if_failed
        )


if __name__ == '__main__':
    luigi.run()
