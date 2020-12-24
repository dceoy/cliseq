#!/usr/bin/env python

import re
import sys
from pathlib import Path

import luigi
from ftarc.task.downloader import DownloadAndProcessResourceFiles
from ftarc.task.picard import CreateSequenceDictionary
from ftarc.task.resource import (FetchReferenceFasta, FetchResourceFile,
                                 FetchResourceVcf)
from luigi.util import requires

from .callcopyratiosegments import PreprocessIntervals
from .core import VclineTask
from .delly import CreateExclusionIntervalListBed
from .msisensorpro import (ScanMicrosatellites,
                           UncompressEvaluationIntervalListBed)


class WritePassingAfOnlyVcf(VclineTask):
    src_path = luigi.Parameter(default='')
    src_url = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    wget = luigi.Parameter(default='wget')
    bgzip = luigi.Parameter(default='bgzip')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(Path(self.src_path or self.src_url).stem).stem
                + '.af-only.vcf.gz'
            )
        )

    def run(self):
        assert bool(self.src_path or self.src_url)
        output_vcf = Path(self.output().path)
        run_id = Path(Path(output_vcf.stem).stem).stem
        message = (
            'Write a passing AF-only VCF' if self.src_path
            else 'Download a VCF file and extract passing AF-only records'
        )
        self.print_log(f'{message}:\t{run_id}')
        dest_dir = output_vcf.parent
        pyscript = Path(__file__).resolve().parent.parent.joinpath(
            'script/extract_af_only_vcf.py'
        )
        self.setup_shell(
            run_id=run_id,
            commands=[
                *(list() if self.src_path else [self.wget]), self.bgzip,
                sys.executable
            ],
            cwd=dest_dir, **self.sh_config
        )
        if self.src_path:
            src_vcf = Path(self.src_path).resolve()
        else:
            src_vcf = dest_dir.joinpath(Path(self.src_url).name)
            self.run_shell(
                args=f'set -e && {self.wget} -qSL {self.src_url} -O {src_vcf}',
                output_files_or_dirs=src_vcf
            )
        self.run_shell(
            args=(
                f'set -e && {self.bgzip}'
                + f' -@ {self.n_cpu} -dc {src_vcf}'
                + f' | {sys.executable} {pyscript} -'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {output_vcf}'
            ),
            input_files_or_dirs=src_vcf, output_files_or_dirs=output_vcf
        )


class PreprocessResources(luigi.Task):
    src_url_dict = luigi.DictParameter()
    dest_dir_path = luigi.Parameter(default='.')
    wget = luigi.Parameter(default='wget')
    bgzip = luigi.Parameter(default='bgzip')
    pbzip2 = luigi.Parameter(default='pbzip2')
    pigz = luigi.Parameter(default='pigz')
    bwa = luigi.Parameter(default='bwa')
    samtools = luigi.Parameter(default='samtools')
    tabix = luigi.Parameter(default='tabix')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    use_bwa_mem2 = luigi.BoolParameter(default=False)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def requires(self):
        return [
            DownloadAndProcessResourceFiles(
                src_urls=[
                    v for k, v in self.src_url_dict.items()
                    if k != 'gnomad_vcf'
                ],
                dest_dir_path=self.dest_dir_path, wget=self.wget,
                bgzip=self.bgzip, pbzip2=self.pbzip2, pigz=self.pigz,
                bwa=self.bwa, samtools=self.samtools, tabix=self.tabix,
                gatk=self.gatk, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                use_bwa_mem2=self.use_bwa_mem2, sh_config=self.sh_config
            ),
            WritePassingAfOnlyVcf(
                src_url=self.src_url_dict['gnomad_vcf'],
                dest_dir_path=self.dest_dir_path, wget=self.wget,
                bgzip=self.bgzip, n_cpu=self.n_cpu, sh_config=self.sh_config
            )
        ]

    def output(self):
        path_dict = self._fetch_input_path_dict()
        fa = Path(path_dict['ref_fa'])
        interval = Path(path_dict['evaluation_interval'])
        gnomad_vcf = Path(path_dict['gnomad_vcf'])
        cnv_blacklist = Path(path_dict['cnv_blacklist'])
        return [
            *self.input(),
            *[
                luigi.LocalTarget(
                    interval.parent.joinpath(interval.stem + s)
                ) for s in [
                    '.bed', '.bed.gz', '.bed.gz.tbi', '.exclusion.bed.gz',
                    '.exclusion.bed.gz.tbi', '.preprocessed.wes.interval_list',
                    '.preprocessed.wgs.interval_list'
                ]
            ],
            *[
                luigi.LocalTarget(
                    gnomad_vcf.parent.joinpath(
                        Path(gnomad_vcf.stem).stem + '.biallelic_snp' + s
                    )
                ) for s in ['.vcf.gz', '.vcf.gz.tbi']
            ],
            *[
                luigi.LocalTarget(
                    cnv_blacklist.parent.joinpath(cnv_blacklist.stem + s)
                ) for s in ['.bed.gz', '.bed.gz.tbi']
            ],
            luigi.LocalTarget(
                fa.parent.joinpath(fa.stem + '.microsatellites.tsv')
            )
        ]

    def run(self):
        path_dict = self._fetch_input_path_dict()
        cf = {
            'pigz': self.pigz, 'pbzip2': self.pbzip2, 'bgzip': self.bgzip,
            'bwa': self.bwa, 'samtools': self.samtools, 'tabix': self.tabix,
            'gatk': self.gatk, 'use_bwa_mem2': self.use_bwa_mem2
        }
        yield [
            CreateExclusionIntervalListBed(
                evaluation_interval_path=path_dict['evaluation_interval'],
                cf=cf, n_cpu=self.n_cpu, sh_config=self.sh_config
            ),
            CreateGnomadBiallelicSnpVCF(
                gnomad_vcf_path=path_dict['gnomad_vcf'],
                ref_fa_path=path_dict['ref_fa'],
                evaluation_interval_path=path_dict['evaluation_interval'],
                cf=cf, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            ),
            CreateCnvBlackListBed(
                cnv_blacklist_path=path_dict['cnv_blacklist'], cf=cf,
                n_cpu=self.n_cpu, sh_config=self.sh_config
            ),
            *[
                PreprocessIntervals(
                    ref_fa_path=path_dict['ref_fa'],
                    evaluation_interval_path=path_dict['evaluation_interval'],
                    cnv_blacklist_path=path_dict['cnv_blacklist'],
                    cf={'exome': bool(i), **cf}, n_cpu=self.n_cpu,
                    memory_mb=self.memory_mb, sh_config=self.sh_config
                ) for i in range(2)
            ],
            ScanMicrosatellites(
                ref_fa_path=path_dict['ref_fa'], cf=cf, n_cpu=self.n_cpu,
                sh_config=self.sh_config
            ),
            UncompressEvaluationIntervalListBed(
                evaluation_interval_path=path_dict['evaluation_interval'],
                cf=cf, n_cpu=self.n_cpu, sh_config=self.sh_config
            )
        ]

    def _fetch_input_path_dict(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        return {
            **{
                k: re.sub(
                    r'\.(gz|bz2)$', '',
                    str(dest_dir.joinpath(Path(self.src_url_dict[k]).name))
                ) for k in ['ref_fa', 'evaluation_interval', 'cnv_blacklist']
            },
            'gnomad_vcf': self.input()[1].path
        }


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
class CreateGnomadBiallelicSnpVCF(VclineTask):
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
                    Path(input_vcf.stem).stem + '.biallelic_snp.vcf.gz' + s
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        input_vcf = Path(self.input()[0][0].path)
        run_id = input_vcf.stem
        self.print_log(f'Create a common biallelic SNP VCF:\t{run_id}')
        gatk = self.cf['gatk']
        tabix = self.cf['tabix']
        fa = Path(self.input()[1][0].path)
        evaluation_interval = Path(self.input()[2].path)
        biallelic_snp_vcf = Path(self.output()[0].path)
        self.setup_shell(
            run_id=run_id, commands=[gatk, tabix], cwd=input_vcf.parent,
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
            output_files_or_dirs=biallelic_snp_vcf
        )
        self.tabix_tbi(tsv_path=biallelic_snp_vcf, tabix=tabix, preset='vcf')


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
                blacklist.parent.joinpath(blacklist.stem + '.bed.gz' + s)
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
