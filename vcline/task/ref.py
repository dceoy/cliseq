#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .base import ShellTask


class FetchReferenceFASTA(ShellTask):
    ref_fa_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            self.ref_fa_paths[0] if (
                len(self.ref_fa_paths) == 1
                and (str(Path(self.ref_fa_paths[0]).parent)
                     == self.cf['ref_dir_path'])
                and self.ref_fa_paths[0].lower().endswith(('.fa', '.fasta'))
            ) else str(
                Path(self.cf['ref_dir_path']).joinpath(
                    '.'.join([
                        re.sub(r'\.(fa|fasta)$', '', Path(p).stem, flags=re.I)
                        for p in self.ref_fa_paths
                    ]) + '.fa'
                )
            )
        )

    def run(self):
        fa_path = self.output().path
        run_id = Path(fa_path).stem
        print_log(f'Create a reference FASTA:\t{run_id}')
        cat = self.cf['cat']
        pigz = self.cf['pigz']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        args = [
            f'{cat} --version',
            f'{pigz} --version',
            f'{pbzip2} --version'
        ]
        for i, p in enumerate(self.ref_fa_paths):
            r = '>' if i == 0 else '>>'
            if p.endswith('.gz'):
                a = f'set -e && {pigz} -p {n_cpu} -dc {p} {r} {fa_path}'
            elif p.endswith('.bz2'):
                a = f'set -e && {pbzip2} -p{n_cpu} -dc {p} {r} {fa_path}'
            else:
                a = 'set -e && {cat} {p} {r} {fa_path}'
            args.append(a)
        self.run_bash(
            args=args, input_files=self.ref_fa_paths, output_files=fa_path
        )


class FetchKnownSiteVCFs(ShellTask):
    known_site_vcf_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 7

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['ref_dir_path']).joinpath(
                        Path(p).name + ('' if p.endswith('.gz') else '.gz')
                    )
                )
            ) for p in self.known_site_vcf_paths
        ]

    def run(self):
        run_id = 'known_site_vcf'
        print_log(f'Create a reference FASTA:\t{run_id}')
        bgzip = self.cf['bgzip']
        n_cpu = self.cf['n_cpu_per_worker']
        vcf_gz_paths = [o.path for o in self.output()]
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        self.run_bash(
            args=[
                f'{bgzip} --version',
                *[
                    'set -e && ' + (
                        f'cp {s} {d}' if s.endswith('.gz') else
                        f'{bgzip} -@ {n_cpu} -c {s} > {d}'
                    ) for s, d in zip(self.known_site_vcf_paths, vcf_gz_paths)
                ]
            ],
            input_files=self.known_site_vcf_paths, output_files=vcf_gz_paths
        )


@requires(FetchKnownSiteVCFs)
class CreateTabixIndices(ShellTask):
    cf = luigi.DictParameter()
    priority = 5

    def output(self):
        return [luigi.LocalTarget(i.path + '.tbi') for i in self.input()]

    def run(self):
        run_id = 'known_site_vcf'
        print_log(f'Create VCF tabix indices:\t{run_id}')
        vcf_gz_paths = [i.path for i in self.input()]
        tabix = self.cf['tabix']
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        self.run_bash(
            args=[
                f'{tabix} --version',
                *[f'set -e && {tabix} -p vcf {p}' for p in vcf_gz_paths]
            ],
            input_files=vcf_gz_paths,
            output_files=[o.path for o in self.output()]
        )


@requires(FetchReferenceFASTA)
class CreateFASTAIndex(ShellTask):
    cf = luigi.DictParameter()
    priority = 8

    def output(self):
        return luigi.LocalTarget(self.input().path + '.fai')

    def run(self):
        fa_path = self.input().path
        run_id = Path(fa_path).stem
        print_log(f'Create a FASTA index:\t{run_id}')
        samtools = self.cf['samtools']
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        self.run_bash(
            args=[
                f'{samtools} 2>&1 | grep -e "Version:"',
                f'set -e && {samtools} faidx {fa_path}'
            ],
            input_files=fa_path, output_files=self.output().path
        )


@requires(FetchReferenceFASTA)
class CreateBWAIndices(ShellTask):
    cf = luigi.DictParameter()
    priority = 9

    def output(self):
        return [
            luigi.LocalTarget(self.input().path + s)
            for s in ['.pac', '.bwt', '.ann', '.amb', '.sa']
        ]

    def run(self):
        fa_path = self.input().path
        run_id = Path(fa_path).stem
        print_log(f'Create BWA indices:\t{run_id}')
        bwa = self.cf['bwa']
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        self.run_bash(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'set -e && {bwa} index {fa_path}'
            ],
            input_files=fa_path,
            output_files=[o.path for o in self.output()]
        )


@requires(FetchReferenceFASTA)
class CreateSequenceDictionary(ShellTask):
    cf = luigi.DictParameter()
    priority = 6

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['ref_dir_path']).joinpath(
                    Path(self.input().path).stem + '.dict'
                )
            )
        )

    def run(self):
        fa_path = self.input().path
        run_id = Path(fa_path).stem
        print_log(f'Create a sequence dictionary:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        self.run_bash(
            args=[
                f'{gatk} --version',
                (
                    'set -e && '
                    + f'{gatk}{gatk_opts} CreateSequenceDictionary'
                    + f' --REFERENCE {fa_path}'
                )
            ],
            input_files=fa_path, output_files=self.output().path
        )


if __name__ == '__main__':
    luigi.run()
