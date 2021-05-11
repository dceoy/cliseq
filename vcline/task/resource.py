#!/usr/bin/env python

import sys
from pathlib import Path

import luigi
from ftarc.task.picard import CreateSequenceDictionary
from ftarc.task.resource import (FetchReferenceFasta, FetchResourceFile,
                                 FetchResourceVcf)
from luigi.util import requires

from .core import VclineTask


class FetchDbsnpVcf(luigi.WrapperTask):
    dbsnp_vcf_path = luigi.Parameter()
    sh_config = luigi.DictParameter(default=dict())
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVcf(
            src_path=self.dbsnp_vcf_path, cf=self.cf, sh_config=self.sh_config
        )

    def output(self):
        return self.input()


class FetchMillsIndelVcf(luigi.WrapperTask):
    mills_indel_vcf_path = luigi.Parameter()
    sh_config = luigi.DictParameter(default=dict())
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVcf(
            src_path=self.mills_indel_vcf_path, cf=self.cf,
            sh_config=self.sh_config
        )

    def output(self):
        return self.input()


class FetchKnownIndelVcf(luigi.WrapperTask):
    known_indel_vcf_path = luigi.Parameter()
    sh_config = luigi.DictParameter(default=dict())
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVcf(
            src_path=self.known_indel_vcf_path, cf=self.cf,
            sh_config=self.sh_config
        )


class FetchEvaluationIntervalList(luigi.WrapperTask):
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 80

    def requires(self):
        return FetchResourceFile(
            src_path=self.evaluation_interval_path, cf=self.cf,
            n_cpu=self.n_cpu, sh_config=self.sh_config
        )

    def output(self):
        return self.input()


@requires(FetchEvaluationIntervalList)
class CreateEvaluationIntervalListBed(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        interval = Path(self.input().path)
        return [
            luigi.LocalTarget(
                interval.parent.joinpath(f'{interval.stem}.bed.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        interval = Path(self.input().path)
        run_id = interval.stem
        self.print_log(f'Create an interval_list BED:\t{run_id}')
        bgzip = self.cf['bgzip']
        tabix = self.cf['tabix']
        bed = Path(self.output()[0].path)
        pyscript = Path(__file__).resolve().parent.parent.joinpath(
            'script/interval_list2bed.py'
        )
        self.setup_shell(
            run_id=run_id, commands=[bgzip, tabix], cwd=interval.parent,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {sys.executable} {pyscript} {interval}'
                + f' | {bgzip} -@ {self.n_cpu} -c > {bed}'
            ),
            input_files_or_dirs=interval, output_files_or_dirs=bed
        )
        self.tabix_tbi(tsv_path=bed, tabix=tabix, preset='bed')


class FetchGnomadVcf(luigi.WrapperTask):
    gnomad_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 80

    def requires(self):
        return FetchResourceVcf(
            src_path=self.gnomad_vcf_path, cf=self.cf, n_cpu=self.n_cpu,
            sh_config=self.sh_config
        )

    def output(self):
        return self.input()


@requires(FetchGnomadVcf, FetchReferenceFasta,
          FetchEvaluationIntervalList, CreateSequenceDictionary)
class CreateGnomadBiallelicSnpVcf(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 90

    def output(self):
        input_vcf = Path(self.input()[0][0].path)
        return [
            luigi.LocalTarget(
                input_vcf.parent.joinpath(
                    Path(input_vcf.stem).stem + f'.biallelic_snp.vcf.gz{s}'
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        input_vcf = Path(self.input()[0][0].path)
        run_id = input_vcf.stem
        self.print_log(f'Create a common biallelic SNP VCF:\t{run_id}')
        gatk = self.cf['gatk']
        fa = Path(self.input()[1][0].path)
        evaluation_interval = Path(self.input()[2].path)
        biallelic_snp_vcf = Path(self.output()[0].path)
        self.setup_shell(
            run_id=run_id, commands=gatk, cwd=input_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {gatk} SelectVariants'
                + f' --variant {input_vcf}'
                + f' --reference {fa}'
                + f' --intervals {evaluation_interval}'
                + f' --output {biallelic_snp_vcf}'
                + ' --select-type-to-include SNP'
                + ' --restrict-alleles-to BIALLELIC'
                + ' --lenient'
            ),
            input_files_or_dirs=[input_vcf, fa, evaluation_interval],
            output_files_or_dirs=[
                biallelic_snp_vcf, f'{biallelic_snp_vcf}.tbi'
            ]
        )


class FetchCnvBlackList(luigi.WrapperTask):
    cnv_blacklist_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 80

    def requires(self):
        return FetchResourceFile(
            src_path=self.cnv_blacklist_path, cf=self.cf, n_cpu=self.n_cpu,
            sh_config=self.sh_config
        )

    def output(self):
        return self.input()


@requires(FetchCnvBlackList)
class CreateCnvBlackListBed(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        blacklist = Path(self.input().path)
        return [
            luigi.LocalTarget(
                blacklist.parent.joinpath(blacklist.stem + f'.bed.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        blacklist = Path(self.input().path)
        run_id = blacklist.stem
        self.print_log(f'Create a blacklist BED:\t{run_id}')
        bgzip = self.cf['bgzip']
        tabix = self.cf['tabix']
        bed = Path(self.output()[0].path)
        self.setup_shell(
            run_id=run_id, commands=[bgzip, tabix], cwd=blacklist.parent,
            **self.sh_config
        )
        pycmd = (
            'from fileinput import input;'
            '[(lambda a, b, c: print(f"{a}\t{int(b) - 1}\t{c}"))'
            '(*s.strip().replace(":", "-").split("-")) for s in input()];'
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {sys.executable}'
                + f' -c \'{pycmd}\' {blacklist}'
                + f' | {bgzip} -@ {self.n_cpu} -c > {bed}'
            ),
            input_files_or_dirs=blacklist, output_files_or_dirs=bed
        )
        self.tabix_tbi(tsv_path=bed, tabix=tabix, preset='bed')


class FetchHapmapVcf(luigi.WrapperTask):
    hapmap_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def requires(self):
        return FetchResourceVcf(
            src_path=self.hapmap_vcf_path, cf=self.cf, n_cpu=self.n_cpu,
            sh_config=self.sh_config
        )

    def output(self):
        return self.input()


class Fetch1000gSnpsVcf(luigi.WrapperTask):
    kg_snps_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def requires(self):
        return FetchResourceVcf(
            src_path=self.kg_snps_vcf_path, cf=self.cf,
            n_cpu=self.n_cpu, sh_config=self.sh_config
        )

    def output(self):
        return self.input()


class CreateIntervalListWithBed(VclineTask):
    bed_path = luigi.Parameter()
    seq_dict_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        return luigi.LocalTarget(
            dest_dir.joinpath(Path(self.bed_path).stem + '.interval_list')
        )

    def run(self):
        interval_list = Path(self.output().path)
        run_id = interval_list.stem
        self.print_log(f'Create an interval_list file:\t{run_id}')
        bed = Path(self.bed_path).resolve()
        seq_dict = Path(self.seq_dict_path).resolve()
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=interval_list.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} BedToIntervalList'
                + f' --INPUT {bed}'
                + f' --SEQUENCE_DICTIONARY {seq_dict}'
                + f' --OUTPUT {interval_list}'
            ),
            input_files_or_dirs=[bed, seq_dict],
            output_files_or_dirs=interval_list
        )


if __name__ == '__main__':
    luigi.run()
