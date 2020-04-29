#!/usr/bin/env python

import re
from configparser import RawConfigParser
from pathlib import Path

import luigi

from .base import BaseTask, ShellTask
from .bcftools import NormalizeVCF, SortVCF
from .ref import (CreateSequenceDictionary, FetchReferenceFASTA,
                  FetchResourceFile)


class AnnotateVariantsWithSnpEff(BaseTask):
    input_vcf_path = luigi.Parameter()
    snpeff_config_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = 10

    def requires(self):
        return [
            FetchResourceFile(
                resource_file_path=self.snpeff_config_path, cf=self.cf
            ),
            FetchReferenceFASTA(ref_fa_path=self.ref_fa_path, cf=self.cf),
            CreateSequenceDictionary(
                ref_fa_path=self.ref_fa_path, cf=self.cf
            )
        ]

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['postproc_snpeff_dir_path']).joinpath(
                        re.sub(
                            r'\.vcf$',
                            '{0}.snpeff.vcf.gz{1}'.format(
                                ('.norm' if self.normalize_vcf else ''), s
                            ),
                            Path(self.input_vcf_path).stem
                        )
                    )
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        yield SnpEff(
            input_vcf_path=self.input_vcf_path,
            fa_path=self.input()[1][0].path,
            snpeff_config_path=self.input()[0].path,
            ref_version=self.cf['ref_version'],
            dest_dir_path=self.cf['postproc_snpeff_dir_path'],
            normalize_vcf=self.normalize_vcf,
            norm_dir_path=self.cf['postproc_bcftools_dir_path'],
            snpeff=self.cf['snpEff'], bcftools=self.cf['bcftools'],
            n_cpu=self.cf['n_cpu_per_worker'],
            memory_mb=self.cf['memory_mb_per_worker'],
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )


class SnpEff(ShellTask):
    input_vcf_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    snpeff_config_path = luigi.Parameter()
    ref_version = luigi.Parameter(default='hg38')
    dest_dir_path = luigi.Parameter(default='.')
    normalize_vcf = luigi.BoolParameter(default=False)
    norm_dir_path = luigi.Parameter(default='')
    snpeff = luigi.Parameter()
    bcftools = luigi.Parameter()
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
                str(
                    Path(self.dest_dir_path).joinpath(
                        re.sub(
                            r'\.vcf$', f'.snpeff.vcf.gz{s}',
                            Path(
                                self.input()[0].path if self.normalize_vcf
                                else self.input_vcf_path
                            ).stem
                        )
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
        tmp_vcf_path = re.sub(r'\.gz$', '', output_vcf_path)
        config = RawConfigParser()
        config.read(self.snpeff_config_path)
        genome_version = [
            o.name for o in Path(config.get('data.dir')).iterdir() if (
                o.name.startswith(
                    {'hg38': 'GRCh38', 'hg19': 'GRCh38'}[self.ref_version]
                ) and o.is_dir()
            )
        ][0]
        self.setup_shell(
            run_id=run_id, log_dir_path=(self.log_dir_path or None),
            commands=[self.snpeff, self.bgzip], cwd=self.dest_dir_path,
            remove_if_failed=self.remove_if_failed,
            quiet=bool(self.log_dir_path)
        )
        self.run_shell(
            args=(
                f'set -e && {self.snpeff} {genome_version} {input_vcf_path}'
                + f' > {tmp_vcf_path}'
            ),
            input_files_or_dirs=input_vcf_path,
            output_files_or_dirs=tmp_vcf_path
        )
        yield SortVCF(
            input_vcf_path=tmp_vcf_path, output_vcf_path=output_vcf_path,
            bcftools=self.bcftools, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            index_vcf=True, remove_input=True, log_dir_path=self.log_dir_path,
            remove_if_failed=self.remove_if_failed,
        )


if __name__ == '__main__':
    luigi.run()
