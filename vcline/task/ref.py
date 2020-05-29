#!/usr/bin/env python

import re
import sys
from itertools import product
from pathlib import Path

import luigi
from luigi.util import requires

from .base import ShellTask
from .samtools import samtools_faidx


class FetchReferenceFASTA(luigi.WrapperTask):
    ref_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        return FetchResourceFASTA(src_path=self.ref_fa_path, cf=self.cf)

    def output(self):
        return self.input()


class FetchResourceFile(ShellTask):
    src_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath(
                re.sub(r'\.(gz|bz2)$', '', Path(self.src_path).name)
            )
        )

    def run(self):
        dest_file = Path(self.output().path)
        run_id = dest_file.stem
        self.print_log(f'Create a resource:\t{run_id}')
        pigz = self.cf['pigz']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[pigz, pbzip2], cwd=dest_file.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        if self.src_path.endswith('.gz'):
            a = f'{pigz} -p {n_cpu} -dc {self.src_path} > {dest_file}'
        elif self.src_path.endswith('.bz2'):
            a = f'{pbzip2} -p{n_cpu} -dc {self.src_path} > {dest_file}'
        else:
            a = f'cp {self.src_path} {dest_file}'
        self.run_shell(
            args=a, input_files_or_dirs=self.src_path,
            output_files_or_dirs=dest_file
        )


@requires(FetchResourceFile)
class FetchResourceFASTA(ShellTask):
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        fa_path = self.input().path
        return [luigi.LocalTarget(fa_path + s) for s in ['', '.fai']]

    def run(self):
        fa = Path(self.input().path)
        run_id = fa.stem
        self.print_log(f'Index FASTA:\t{run_id}')
        samtools = self.cf['samtools']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=samtools, cwd=fa.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        samtools_faidx(shelltask=self, samtools=samtools, fa_path=fa)


@requires(FetchReferenceFASTA)
class CreateSequenceDictionary(ShellTask):
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        fa = Path(self.input()[0].path)
        return luigi.LocalTarget(fa.parent.joinpath(f'{fa.stem}.dict'))

    def run(self):
        fa = Path(self.input()[0].path)
        run_id = fa.stem
        self.print_log(f'Create a sequence dictionary:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        seq_dict_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=fa.parent, remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} CreateSequenceDictionary'
                + f' --REFERENCE {fa}'
                + f' --OUTPUT {seq_dict_path}'
            ),
            input_files_or_dirs=fa, output_files_or_dirs=seq_dict_path
        )


class FetchResourceVCF(ShellTask):
    src_path = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        dest_vcf = Path(self.cf['ref_dir_path']).joinpath(
            re.sub(r'\.(gz|bgz)$', '.gz', Path(self.src_path).name)
        )
        return [luigi.LocalTarget(f'{dest_vcf}{s}') for s in ['', '.tbi']]

    def run(self):
        dest_vcf = Path(self.output()[0].path)
        run_id = Path(dest_vcf.stem).stem
        self.print_log(f'Create a VCF:\t{run_id}')
        bgzip = self.cf['bgzip']
        tabix = self.cf['tabix']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[bgzip, tabix], cwd=dest_vcf.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                'set -e && ' + (
                    f'cp {self.src_path} {dest_vcf}'
                    if self.src_path.endswith(('.gz', '.bgz')) else
                    f'{bgzip} -@ {n_cpu} -c {self.src_path} > {dest_vcf}'
                )
            ),
            input_files_or_dirs=self.src_path, output_files_or_dirs=dest_vcf
        )
        _tabix(shelltask=self, tabix=tabix, tsv_path=dest_vcf, preset='vcf')


def _tabix(shelltask, tabix, tsv_path, preset='bed'):
    shelltask.run_shell(
        args=f'set -e && {tabix} --preset {preset} {tsv_path}',
        input_files_or_dirs=tsv_path, output_files_or_dirs=f'{tsv_path}.tbi'
    )


class FetchDbsnpVCF(luigi.WrapperTask):
    dbsnp_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVCF(src_path=self.dbsnp_vcf_path, cf=self.cf)

    def output(self):
        return self.input()


class FetchMillsIndelVCF(luigi.WrapperTask):
    mills_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVCF(src_path=self.mills_indel_vcf_path, cf=self.cf)

    def output(self):
        return self.input()


class FetchKnownIndelVCF(luigi.WrapperTask):
    known_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVCF(src_path=self.known_indel_vcf_path, cf=self.cf)

    def output(self):
        return self.input()


@requires(FetchReferenceFASTA)
class CreateBWAIndices(ShellTask):
    cf = luigi.DictParameter()
    priority = 100

    def output(self):
        fa_path = self.input()[0].path
        return [
            luigi.LocalTarget(f'{fa_path}.{s}')
            for s in ['pac', 'bwt', 'ann', 'amb', 'sa']
        ]

    def run(self):
        fa = Path(self.input()[0].path)
        run_id = fa.stem
        self.print_log(f'Create BWA indices:\t{run_id}')
        bwa = self.cf['bwa']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bwa, cwd=fa.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=f'set -e && {bwa} index {fa}', input_files_or_dirs=fa,
            output_files_or_dirs=[o.path for o in self.output()]
        )


class FetchEvaluationIntervalList(luigi.WrapperTask):
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 80

    def requires(self):
        return FetchResourceFile(
            src_path=self.evaluation_interval_path, cf=self.cf
        )

    def output(self):
        return self.input()


@requires(FetchEvaluationIntervalList)
class CreateEvaluationIntervalListBED(ShellTask):
    cf = luigi.DictParameter()
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
        n_cpu = self.cf['n_cpu_per_worker']
        bed_path = self.output()[0].path
        pyscript_path = str(
            Path(__file__).parent.parent.joinpath(
                'script/interval_list2bed.py'
            ).resolve()
        )
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[bgzip, tabix], cwd=interval.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {sys.executable} {pyscript_path}'
                + f' {interval}'
                + f' | {bgzip} -@ {n_cpu} -c > {bed_path}'
            ),
            input_files_or_dirs=interval, output_files_or_dirs=bed_path
        )
        _tabix(shelltask=self, tabix=tabix, tsv_path=bed_path, preset='bed')


@requires(CreateEvaluationIntervalListBED, FetchReferenceFASTA)
class CreateExclusionIntervalListBED(ShellTask):
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return [
            luigi.LocalTarget(p + s) for p, s in product(
                [
                    re.sub(
                        r'\.bed\.gz$', '.exclusion.bed.gz',
                        self.input()[0][0].path
                    ),
                    (self.input()[1][0].path + '.bed.gz')
                ],
                ['', '.tbi']
            )
        ]

    def run(self):
        evaluation_bed = Path(self.input()[0][0].path)
        run_id = Path(evaluation_bed.stem).stem
        self.print_log(f'Create an exclusion interval_list BED:\t{run_id}')
        bedtools = self.cf['bedtools']
        bgzip = self.cf['bgzip']
        tabix = self.cf['tabix']
        fai_path = self.input()[1][1].path
        n_cpu = self.cf['n_cpu_per_worker']
        exclusion_bed_path = self.output()[0].path
        genome_bed_path = self.output()[2].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[bgzip, tabix], cwd=evaluation_bed.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {sys.executable}'
                + ' -c \'{}\''.format(
                    'from fileinput import input; '
                    '[print("{0}\\t0\\t{1}".format(*s.split()[:2]))'
                    ' for s in input()];'
                ) + f' {fai_path}'
                + f' | {bgzip} -@ {n_cpu} -c > {genome_bed_path}'
            ),
            input_files_or_dirs=fai_path, output_files_or_dirs=genome_bed_path
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {bedtools} subtract'
                + f' -a {genome_bed_path} -b {evaluation_bed}'
                + f' | {bgzip} -@ {n_cpu} -c > {exclusion_bed_path}'
            ),
            input_files_or_dirs=[genome_bed_path, evaluation_bed],
            output_files_or_dirs=exclusion_bed_path
        )
        for p in [genome_bed_path, exclusion_bed_path]:
            _tabix(shelltask=self, tabix=tabix, tsv_path=p, preset='bed')


class FetchHapmapVCF(luigi.WrapperTask):
    hapmap_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVCF(src_path=self.hapmap_vcf_path, cf=self.cf)

    def output(self):
        return self.input()


class FetchGnomadVCF(luigi.WrapperTask):
    gnomad_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 80

    def requires(self):
        return FetchResourceVCF(src_path=self.gnomad_vcf_path, cf=self.cf)

    def output(self):
        return self.input()


@requires(FetchGnomadVCF, FetchReferenceFASTA,
          FetchEvaluationIntervalList, CreateSequenceDictionary)
class CreateGnomadBiallelicSnpVCF(ShellTask):
    cf = luigi.DictParameter()
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
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        tabix = self.cf['tabix']
        fa_path = self.input()[1][0].path
        evaluation_interval_path = self.input()[2].path
        biallelic_snp_vcf_path = self.output()[0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, tabix], cwd=input_vcf.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} SelectVariants'
                + f' --variant {input_vcf}'
                + f' --reference {fa_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {biallelic_snp_vcf_path}'
                + ' --select-type-to-include SNP'
                + ' --restrict-alleles-to BIALLELIC'
                + ' --lenient'
            ),
            input_files_or_dirs=[input_vcf, fa_path, evaluation_interval_path],
            output_files_or_dirs=biallelic_snp_vcf_path
        )
        _tabix(
            shelltask=self, tabix=tabix, tsv_path=biallelic_snp_vcf_path,
            preset='vcf'
        )


class FetchCnvBlackList(luigi.WrapperTask):
    cnv_blacklist_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 80

    def requires(self):
        return FetchResourceFile(
            src_path=self.cnv_blacklist_path, cf=self.cf
        )

    def output(self):
        return self.input()


@requires(FetchCnvBlackList)
class CreateCnvBlackListBED(ShellTask):
    cf = luigi.DictParameter()
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
        n_cpu = self.cf['n_cpu_per_worker']
        bed_path = self.output()[0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[bgzip, tabix], cwd=blacklist.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {sys.executable}'
                + ' -c \'{}\''.format(
                    'from fileinput import input;'
                    '[(lambda a, b, c: print(f"{a}\t{int(b) - 1}\t{c}"))'
                    '(*s.strip().replace(":", "-").split("-"))'
                    ' for s in input()];'
                ) + f' {blacklist}'
                + f' | {bgzip} -@ {n_cpu} -c > {bed_path}'
            ),
            input_files_or_dirs=blacklist, output_files_or_dirs=bed_path
        )
        _tabix(shelltask=self, tabix=tabix, tsv_path=bed_path, preset='bed')


@requires(FetchEvaluationIntervalList, FetchCnvBlackList,
          FetchReferenceFASTA, CreateSequenceDictionary)
class PreprocessIntervals(ShellTask):
    cf = luigi.DictParameter()
    wes_param = luigi.DictParameter(default={'bin-length': 0, 'padding': 250})
    wgs_param = luigi.DictParameter(default={'bin-length': 1000, 'padding': 0})
    priority = 60

    def output(self):
        evaluation_interval = Path(self.input()[0].path)
        return luigi.LocalTarget(
            evaluation_interval.parent.joinpath(
                '{0}.preprocessed.w{1}s.interval_list'.format(
                    Path(self.input()[0].path).stem,
                    ('e' if self.cf['exome'] else 'g')
                )
            )
        )

    def run(self):
        evaluation_interval = Path(self.input()[0].path)
        run_id = evaluation_interval
        self.print_log(f'Prepares bins for coverage collection:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        cnv_blacklist_path = self.input()[1].path
        fa_path = self.input()[2][0].path
        seq_dict_path = self.input()[3].path
        preprocessed_interval_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=evaluation_interval.parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} PreprocessIntervals'
                + f' --intervals {evaluation_interval}'
                + f' --exclude-intervals {cnv_blacklist_path}'
                + f' --sequence-dictionary {seq_dict_path}'
                + f' --reference {fa_path}'
                + ''.join([
                    f' --{k} {v}' for k, v in (
                        self.wes_param if self.cf['exome'] else self.wgs_param
                    ).items()
                ])
                + ' --interval-merging-rule OVERLAPPING_ONLY'
                + f' --output {preprocessed_interval_path}'
            ),
            input_files_or_dirs=[
                evaluation_interval, cnv_blacklist_path, seq_dict_path, fa_path
            ],
            output_files_or_dirs=preprocessed_interval_path
        )


class UncompressBgzipFiles(ShellTask):
    bgz_paths = luigi.ListParameter()
    dest_dir_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        dest_dir = Path(self.dest_dir_path)
        return [
            luigi.LocalTarget(dest_dir.joinpath(Path(p).stem))
            for p in self.bgz_paths
        ]

    def run(self):
        run_id = Path(self.dest_dir_path).name
        self.print_log(f'Uncompress bgzip files:\t{run_id}')
        bgzip = self.cf['bgzip']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bgzip, cwd=self.dest_dir_path,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        for p, o in zip(self.bgz_paths, self.output()):
            self.run_shell(
                args=f'set -e && {bgzip} -@ {n_cpu} -dc {p} > {o.path}',
                input_files_or_dirs=p, output_files_or_dirs=o.path
            )


class CreateSymlinks(ShellTask):
    src_paths = luigi.ListParameter()
    dest_dir_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        dest_dir = Path(self.dest_dir_path)
        return [
            luigi.LocalTarget(dest_dir.joinpath(Path(p).name))
            for p in self.src_paths
        ]

    def run(self):
        run_id = Path(self.dest_dir_path).name
        self.print_log(f'Create a symlink:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=self.dest_dir_path,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        for p, o in zip(self.src_paths, self.output()):
            self.run_shell(
                args=f'ln -s {p} {o.path}',
                input_files_or_dirs=p, output_files_or_dirs=o.path
            )


if __name__ == '__main__':
    luigi.run()
