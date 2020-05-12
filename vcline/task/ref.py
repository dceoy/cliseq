#!/usr/bin/env python

import re
import sys
from pathlib import Path

import luigi
from luigi.util import requires

from .base import ShellTask
from .samtools import SamtoolsFaidx, Tabix


class FetchReferenceFASTA(luigi.WrapperTask):
    ref_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        return FetchResourceFASTA(
            resource_file_path=self.ref_fa_path, cf=self.cf
        )

    def output(self):
        return self.input()


class FetchResourceFile(ShellTask):
    resource_file_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath(
                re.sub(r'\.(gz|bz2)$', '', Path(self.resource_file_path).name)
            )
        )

    def run(self):
        dest_path = self.output().path
        run_id = Path(dest_path).stem
        self.print_log(f'Create a resource:\t{run_id}')
        src_path = self.resource_file_path
        pigz = self.cf['pigz']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[pigz, pbzip2], cwd=self.cf['ref_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        if src_path.endswith('.gz'):
            a = f'{pigz} -p {n_cpu} -dc {src_path} > {dest_path}'
        elif src_path.endswith('.bz2'):
            a = f'{pbzip2} -p{n_cpu} -dc {src_path} > {dest_path}'
        else:
            a = f'cp {src_path} {dest_path}'
        self.run_shell(
            args=a, input_files_or_dirs=src_path,
            output_files_or_dirs=dest_path
        )


@requires(FetchReferenceFASTA)
class CreateSequenceDictionary(ShellTask):
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath(
                Path(self.input()[0].path).stem + '.dict'
            )
        )

    def run(self):
        fa_path = self.input()[0].path
        run_id = Path(fa_path).stem
        self.print_log(f'Create a sequence dictionary:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        seq_dict_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['ref_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} CreateSequenceDictionary'
                + f' --REFERENCE {fa_path}'
                + f' --OUTPUT {seq_dict_path}'
            ),
            input_files_or_dirs=fa_path, output_files_or_dirs=seq_dict_path
        )


@requires(FetchResourceFile)
class FetchResourceFASTA(luigi.Task):
    resource_file_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        fa_path = self.input().path
        return [luigi.LocalTarget(f'{fa_path}{s}') for s in ['', '.fai']]

    def run(self):
        yield SamtoolsFaidx(
            fa_path=self.input().path, samtools=self.cf['samtools'],
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )


class FetchResourceVCF(ShellTask):
    resource_vcf_path = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['ref_dir_path']).joinpath(
                    re.sub(
                        r'\.(gz|bgz)$', f'.{s}',
                        Path(self.resource_vcf_path).name
                    )
                )
            ) for s in ['gz', 'gz.tbi']
        ]

    def run(self):
        dest_vcf_path = self.output()[0].path
        run_id = Path(Path(dest_vcf_path).stem).stem
        self.print_log(f'Create a VCF:\t{run_id}')
        bgzip = self.cf['bgzip']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bgzip, cwd=self.cf['ref_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
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
            input_files_or_dirs=self.resource_vcf_path,
            output_files_or_dirs=dest_vcf_path
        )
        yield Tabix(
            tsv_path=dest_vcf_path, tabix=self.cf['tabix'], preset='vcf',
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )


class FetchDbsnpVCF(luigi.WrapperTask):
    dbsnp_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVCF(
            resource_vcf_path=self.dbsnp_vcf_path, cf=self.cf
        )

    def output(self):
        return self.input()


class FetchMillsIndelVCF(luigi.WrapperTask):
    mills_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVCF(
            resource_vcf_path=self.mills_indel_vcf_path, cf=self.cf
        )

    def output(self):
        return self.input()


class FetchKnownIndelVCF(luigi.WrapperTask):
    known_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVCF(
            resource_vcf_path=self.known_indel_vcf_path, cf=self.cf
        )

    def output(self):
        return self.input()


@requires(FetchReferenceFASTA)
class CreateBWAIndices(ShellTask):
    cf = luigi.DictParameter()
    priority = 100

    def output(self):
        return [
            luigi.LocalTarget(self.input()[0].path + s)
            for s in ['.pac', '.bwt', '.ann', '.amb', '.sa']
        ]

    def run(self):
        fa_path = self.input()[0].path
        run_id = Path(fa_path).stem
        self.print_log(f'Create BWA indices:\t{run_id}')
        bwa = self.cf['bwa']
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bwa, cwd=self.cf['ref_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=f'set -e && {bwa} index {fa_path}',
            input_files_or_dirs=fa_path,
            output_files_or_dirs=[o.path for o in self.output()]
        )


class FetchEvaluationIntervalList(luigi.WrapperTask):
    evaluation_interval_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 80

    def requires(self):
        return FetchResourceFile(
            resource_file_path=self.evaluation_interval_path, cf=self.cf
        )

    def output(self):
        return self.input()


@requires(FetchEvaluationIntervalList, FetchReferenceFASTA,
          CreateSequenceDictionary)
class SplitEvaluationIntervals(ShellTask):
    scatter_count = luigi.IntParameter(default=2)
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['ref_dir_path']).joinpath(
                    Path(self.input()[0].path).stem + '.split'
                ).joinpath(f'{i:04d}-scattered.interval_list')
            ) for i in range(self.scatter_count)
        ]

    def run(self):
        input_interval_path = self.input()[0].path
        run_id = Path(input_interval_path).stem
        self.print_log(f'Split an evaluation interval list:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        fa_path = self.input()[1][0].path
        output_interval_paths = [o.path for o in self.output()]
        output_dir = Path(output_interval_paths[0]).parent
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=output_dir, remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                'set -e && '
                + f'{gatk}{gatk_opts} SplitIntervals'
                + f' --reference {fa_path}'
                + f' --intervals {input_interval_path}'
                + f' --scatter-count {self.scatter_count}'
                + f' --output {output_dir}'
            ),
            input_files_or_dirs=[input_interval_path, fa_path],
            output_files_or_dirs=output_interval_paths
        )


class PrepareEvaluationIntervals(luigi.WrapperTask):
    evaluation_interval_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

    def requires(self):
        return (
            SplitEvaluationIntervals(
                evaluation_interval_path=self.evaluation_interval_path,
                ref_fa_path=self.ref_fa_path,
                scatter_count=self.cf['n_worker'], cf=self.cf
            ) if self.cf['n_worker'] > 1 else [
                FetchEvaluationIntervalList(
                    evaluation_interval_path=self.evaluation_interval_path,
                    cf=self.cf
                )
            ]
        )

    def output(self):
        return self.input()


@requires(FetchEvaluationIntervalList)
class CreateEvaluationIntervalListBED(luigi.Task):
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return [
            luigi.LocalTarget(self.input().path + f'.bed.gz{s}')
            for s in ['', '.tbi']
        ]

    def run(self):
        yield IntervalList2BED(
            interval_list_path=self.input().path, cf=self.cf
        )


class IntervalList2BED(ShellTask):
    interval_list_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return [
            luigi.LocalTarget(self.interval_list_path + f'.bed.gz{s}')
            for s in ['', '.tbi']
        ]

    def run(self):
        run_id = Path(self.interval_list_path).stem
        self.print_log(f'Create an interval_list BED:\t{run_id}')
        bgzip = self.cf['bgzip']
        n_cpu = self.cf['n_cpu_per_worker']
        interval_bed_path = self.output()[0].path
        pyscript_path = str(
            Path(__file__).parent.parent.joinpath(
                'script/interval_list2bed.py'
            ).resolve()
        )
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=self.cf['bgzip'], cwd=Path(interval_bed_path).parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && '
                + f'{sys.executable} {pyscript_path} {self.interval_list_path}'
                + f' | {bgzip} -@ {n_cpu} -c > {interval_bed_path}'
            ),
            input_files_or_dirs=self.interval_list_path,
            output_files_or_dirs=interval_bed_path
        )
        yield Tabix(
            tsv_path=interval_bed_path, tabix=self.cf['tabix'], preset='bed',
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )


@requires(CreateEvaluationIntervalListBED, FetchReferenceFASTA)
class CreateExclusionIntervalListBED(ShellTask):
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return (
            [
                luigi.LocalTarget(
                    re.sub(
                        r'\.bed\.gz$', f'.exclusion.bed.gz{s}',
                        self.input()[0][0].path
                    )
                ) for s in ['', '.tbi']
            ] + [
                luigi.LocalTarget(self.input()[1][0].path + f'.bed.gz{s}')
                for s in ['', '.tbi']
            ]
        )

    def run(self):
        evaluation_bed_path = self.input()[0][0].path
        run_id = Path(Path(evaluation_bed_path).stem).stem
        self.print_log(f'Create an exclusion interval_list BED:\t{run_id}')
        bedtools = self.cf['bedtools']
        bgzip = self.cf['bgzip']
        fai_path = self.input()[1][1].path
        n_cpu = self.cf['n_cpu_per_worker']
        exclusion_bed_path = self.output()[0].path
        genome_bed_path = self.output()[2].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bgzip, cwd=self.cf['ref_dir_path'],
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
                + f' -a {genome_bed_path} -b {evaluation_bed_path}'
                + f' | {bgzip} -@ {n_cpu} -c > {exclusion_bed_path}'
            ),
            input_files_or_dirs=[genome_bed_path, evaluation_bed_path],
            output_files_or_dirs=exclusion_bed_path
        )
        yield [
            Tabix(
                tsv_path=p, tabix=self.cf['tabix'], preset='bed',
                log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed'],
                quiet=self.cf['quiet']
            ) for p in [genome_bed_path, exclusion_bed_path]
        ]


class FetchHapmapVCF(luigi.WrapperTask):
    hapmap_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 70

    def requires(self):
        return FetchResourceVCF(
            resource_vcf_path=self.hapmap_vcf_path, cf=self.cf
        )

    def output(self):
        return self.input()


class FetchGnomadVCF(luigi.WrapperTask):
    gnomad_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 80

    def requires(self):
        return FetchResourceVCF(
            resource_vcf_path=self.gnomad_vcf_path, cf=self.cf
        )

    def output(self):
        return self.input()


@requires(FetchGnomadVCF, FetchReferenceFASTA,
          FetchEvaluationIntervalList, CreateSequenceDictionary)
class CreateGnomadBiallelicSnpVCF(ShellTask):
    cf = luigi.DictParameter()
    priority = 90

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['ref_dir_path']).joinpath(
                    Path(Path(self.input()[0][0].path).stem).stem
                    + f'.biallelic_snp.vcf.{s}'
                )
            ) for s in ['gz', 'gz.tbi']
        ]

    def run(self):
        input_vcf_path = self.input()[0][0].path
        run_id = Path(input_vcf_path).stem
        self.print_log(f'Create a common biallelic SNP VCF:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        fa_path = self.input()[1][0].path
        evaluation_interval_path = self.input()[2].path
        biallelic_snp_vcf_path = self.output()[0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['ref_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} SelectVariants'
                + f' --variant {input_vcf_path}'
                + f' --reference {fa_path}'
                + f' --intervals {evaluation_interval_path}'
                + f' --output {biallelic_snp_vcf_path}'
                + ' --select-type-to-include SNP'
                + ' --restrict-alleles-to BIALLELIC'
                + ' --lenient'
            ),
            input_files_or_dirs=[
                input_vcf_path, fa_path, evaluation_interval_path
            ],
            output_files_or_dirs=biallelic_snp_vcf_path
        )
        yield Tabix(
            tsv_path=biallelic_snp_vcf_path, tabix=self.cf['tabix'],
            preset='vcf', log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )


class ExtractTarFile(ShellTask):
    tar_path = luigi.Parameter()
    cf = luigi.DictParameter()
    recursive = luigi.BoolParameter(default=True)
    priority = 70

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath(
                re.sub(
                    r'\.tar\.(gz|bz2)$', '', Path(self.tar_path).name
                )
            )
        )

    def run(self):
        dest_path = self.output().path
        run_id = Path(dest_path).name
        self.print_log(f'Create a resource:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=self.cf['ref_dir_path']
        )
        self.run_shell(
            args=f'tar -xvf {self.tar_path}',
            input_files_or_dirs=self.tar_path, output_files_or_dirs=dest_path
        )
        if self.recursive and Path(dest_path).is_dir():
            for o in Path(dest_path).iterdir():
                if o.name.endswith(('.tar.gz', '.tar.bz2')):
                    p = str(o)
                    self.run_shell(
                        args=f'tar -xvf {p}',
                        cwd=dest_path, input_files_or_dirs=p,
                        output_files_or_dirs=re.sub(r'\.tar\.(gz|bz2)$', '', p)
                    )


class FetchCnvBlackList(luigi.WrapperTask):
    cnv_black_list_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 80

    def requires(self):
        return FetchResourceFile(
            resource_file_path=self.cnv_black_list_path, cf=self.cf
        )

    def output(self):
        return self.input()


@requires(FetchCnvBlackList)
class CreateCnvBlackListBED(ShellTask):
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return [
            luigi.LocalTarget(self.input().path + f'.bed.gz{s}')
            for s in ['', '.tbi']
        ]

    def run(self):
        blacklist_path = self.input().path
        run_id = Path(blacklist_path).stem
        self.print_log(f'Create a blacklist BED:\t{run_id}')
        bgzip = self.cf['bgzip']
        n_cpu = self.cf['n_cpu_per_worker']
        bed_path = self.output()[0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bgzip, cwd=self.cf['ref_dir_path'],
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
                ) + f' {blacklist_path}'
                + f' | {bgzip} -@ {n_cpu} -c > {bed_path}'
            ),
            input_files_or_dirs=blacklist_path, output_files_or_dirs=bed_path
        )
        yield Tabix(
            tsv_path=bed_path, tabix=self.cf['tabix'], preset='bed',
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )


@requires(FetchEvaluationIntervalList, FetchCnvBlackList,
          FetchReferenceFASTA, CreateSequenceDictionary)
class PreprocessIntervals(ShellTask):
    cf = luigi.DictParameter()
    wes_param = luigi.DictParameter(default={'bin-length': 0, 'padding': 250})
    wgs_param = luigi.DictParameter(default={'bin-length': 1000, 'padding': 0})
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath(
                '{0}.preprocessed.w{1}s.interval_list'.format(
                    Path(self.input()[0].path).stem,
                    ('e' if self.cf['exome'] else 'g')
                )
            )
        )

    def run(self):
        evaluation_interval_path = self.input()[0].path
        run_id = Path(evaluation_interval_path).stem
        self.print_log(f'Prepares bins for coverage collection:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        cnv_black_list_path = self.input()[1].path
        fa_path = self.input()[2][0].path
        seq_dict_path = self.input()[3].path
        preprocessed_interval_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'], commands=gatk,
            cwd=self.cf['ref_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} PreprocessIntervals'
                + f' --intervals {evaluation_interval_path}'
                + f' --exclude-intervals {cnv_black_list_path}'
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
                evaluation_interval_path, cnv_black_list_path,
                seq_dict_path, fa_path
            ],
            output_files_or_dirs=preprocessed_interval_path
        )


class UncompressBgzipFiles(ShellTask):
    bgz_paths = luigi.ListParameter()
    dest_dir_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return [
            luigi.LocalTarget(Path(self.dest_dir_path).joinpath(Path(p).stem))
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
        return [
            luigi.LocalTarget(Path(self.dest_dir_path).joinpath(Path(p).name))
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
