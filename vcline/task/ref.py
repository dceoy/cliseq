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
    priority = 100

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
        pigz = self.cf['pigz']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[pigz, pbzip2], cwd=self.cf['ref_dir_path']
        )
        args = list()
        for i, p in enumerate(self.ref_fa_paths):
            r = '>' if i == 0 else '>>'
            if p.endswith('.gz'):
                a = f'{pigz} -p {n_cpu} -dc {p} {r} {fa_path}'
            elif p.endswith('.bz2'):
                a = f'{pbzip2} -p{n_cpu} -dc {p} {r} {fa_path}'
            else:
                a = 'cat {p} {r} {fa_path}'
            args.append(f'set -e && {a}')
        self.run_shell(
            args=args, input_files=self.ref_fa_paths, output_files=fa_path
        )


class FetchResourceVCF(ShellTask):
    resource_vcf_path = luigi.ListParameter()
    cf = luigi.DictParameter()

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['ref_dir_path']).joinpath(
                        re.sub(
                            r'\.(gz|bgz)$', f'.{s}',
                            Path(self.resource_vcf_path).name
                        )
                    )
                )
            ) for s in ['gz', 'gz.tbi']
        ]

    def run(self):
        dest_vcf_path = self.output()[0].path
        run_id = Path(Path(dest_vcf_path).stem).stem
        print_log(f'Create a VCF:\t{run_id}')
        bgzip = self.cf['bgzip']
        tabix = self.cf['tabix']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[bgzip, tabix], cwd=self.cf['ref_dir_path']
        )
        self.run_shell(
            args=(
                'set -e && ' + (
                    f'cp {self.resource_vcf_path} {dest_vcf_path}'
                    if self.resource_vcf_path.endswith(('.gz', '.bgz')) else (
                        f'{bgzip} -@ {n_cpu} -c'
                        + f' {self.resource_vcf_path} > {dest_vcf_path}'
                    )
                )
            ),
            input_files=self.resource_vcf_path, output_files=dest_vcf_path
        )
        self.run_shell(
            args=f'set -e && {tabix} -p vcf {dest_vcf_path}',
            input_files=dest_vcf_path, output_files=f'{dest_vcf_path}.tbi'
        )


class FetchDbsnpVCF(luigi.WrapperTask):
    dbsnp_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 80

    def requires(self):
        return FetchResourceVCF(
            resource_vcf_path=self.dbsnp_vcf_path, cf=self.cf
        )

    def output(self):
        return self.input()


class FetchKnownIndelVCFs(luigi.WrapperTask):
    known_indel_vcf_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 80

    def requires(self):
        return [
            FetchResourceVCF(resource_vcf_path=p, cf=self.cf)
            for p in self.known_indel_vcf_paths
        ]

    def output(self):
        return self.input()


class FetchResourceFile(ShellTask):
    resource_file_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 80

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['ref_dir_path']).joinpath(
                    re.sub(
                        r'\.(gz|bz2)$', '', Path(self.resource_file_path).name
                    )
                )
            )
        )

    def run(self):
        dest_path = self.output().path
        run_id = Path(dest_path).stem
        print_log(f'Create a resource:\t{run_id}')
        src_path = self.resource_file_path
        pigz = self.cf['pigz']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[pigz, pbzip2], cwd=self.cf['ref_dir_path']
        )
        if src_path.endswith('.gz'):
            a = f'{pigz} -p {n_cpu} -dc {src_path} > {dest_path}'
        elif src_path.endswith('.bz2'):
            a = f'{pbzip2} -p{n_cpu} -dc {src_path} > {dest_path}'
        else:
            a = f'cp {src_path} {dest_path}'
        self.run_shell(args=a, input_files=src_path, output_files=dest_path)


@requires(FetchReferenceFASTA)
class CreateEvaluationIntervalList(ShellTask):
    cf = luigi.DictParameter()
    priority = 90

    def output(self):
        return luigi.LocalTarget(self.input().path + '.interval_list')

    def run(self):
        fa_path = self.input().path
        run_id = Path(fa_path).stem
        print_log(f'Create an evaluation interval list:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        interval_list_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['ref_dir_path']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} ScatterIntervalsByNs'
                + f' --REFERENCE {fa_path}'
                + f' --OUTPUT {interval_list_path}'
                + ' --OUTPUT_TYPE=ACGT'
            ),
            input_files=fa_path, output_files=interval_list_path
        )


@requires(FetchReferenceFASTA)
class CreateFASTAIndex(ShellTask):
    cf = luigi.DictParameter()
    priority = 80

    def output(self):
        return luigi.LocalTarget(self.input().path + '.fai')

    def run(self):
        fa_path = self.input().path
        run_id = Path(fa_path).stem
        print_log(f'Create a FASTA index:\t{run_id}')
        samtools = self.cf['samtools']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=samtools, cwd=self.cf['ref_dir_path']
        )
        self.run_shell(
            args=f'set -e && {samtools} faidx {fa_path}',
            input_files=fa_path, output_files=self.output().path
        )


@requires(FetchReferenceFASTA)
class CreateBWAIndices(ShellTask):
    cf = luigi.DictParameter()
    priority = 100

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
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bwa, cwd=self.cf['ref_dir_path']
        )
        self.run_shell(
            args=f'set -e && {bwa} index {fa_path}',
            input_files=fa_path,
            output_files=[o.path for o in self.output()]
        )


@requires(FetchReferenceFASTA)
class CreateSequenceDictionary(ShellTask):
    cf = luigi.DictParameter()
    priority = 80

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
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['ref_dir_path']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} CreateSequenceDictionary'
                + f' --REFERENCE {fa_path}'
            ),
            input_files=fa_path, output_files=self.output().path
        )


class PrepareGermlineResourceVCFs(luigi.WrapperTask):
    dbsnp_vcf_path = luigi.Parameter()
    hapmap_vcf_path = luigi.Parameter()
    omni_vcf_path = luigi.Parameter()
    snp_1000g_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 40

    def requires(self):
        return {
            'dbsnp': FetchDbsnpVCF(
                dbsnp_vcf_path=self.dbsnp_vcf_path, cf=self.cf
            ),
            'hapmap': FetchResourceVCF(
                resource_vcf_path=self.hapmap_vcf_path, cf=self.cf
            ),
            'omni': FetchResourceVCF(
                resource_vcf_path=self.omni_vcf_path, cf=self.cf
            ),
            'snp_1000g': FetchResourceVCF(
                resource_vcf_path=self.snp_1000g_vcf_path, cf=self.cf
            )
        }

    def output(self):
        return self.input()


class CreateGnomadSelectedVCF(ShellTask):
    gnomad_vcf_path = luigi.Parameter()
    ref_fa_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 50

    def requires(self):
        return [
            FetchResourceVCF(
                resource_vcf_path=self.gnomad_vcf_path, cf=self.cf
            ),
            FetchReferenceFASTA(
                ref_fa_paths=self.ref_fa_paths, cf=self.cf
            ),
            CreateEvaluationIntervalList(
                ref_fa_paths=self.ref_fa_paths, cf=self.cf
            )
        ]

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['call_dir_path']).joinpath(
                        Path(self.input()[0][0].path).stem
                        + f'.selected.vcf.{s}'
                    )
                )
            ) for s in ['gz', 'gz.tbi']
        ]

    def run(self):
        gnomad_vcf_path = self.input()[0][0].path
        run_id = '.'.join(Path(gnomad_vcf_path).name.split('.')[:-2])
        print_log(f'Create a gnomAD VCF:\t{run_id}')
        tabix = self.cf['tabix']
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        fa_path = self.input()[1].path
        evaluation_interval_path = self.input()[2].path
        selected_vcf_path = self.output()[0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[tabix, gatk], cwd=self.cf['ref_dir_path']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} SelectVariants'
                + f' --reference {fa_path}'
                + f' --variant {gnomad_vcf_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {selected_vcf_path}'
                + ' --lenient'
            ),
            input_files=[gnomad_vcf_path, fa_path, evaluation_interval_path],
            output_files=selected_vcf_path
        )
        self.run_shell(
            args=f'set -e && {tabix} -p vcf {selected_vcf_path}',
            input_files=selected_vcf_path,
            output_files=f'{selected_vcf_path}.tbi'
        )


if __name__ == '__main__':
    luigi.run()
