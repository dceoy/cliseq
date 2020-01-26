#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import parse_fq_id, print_log
from .base import ShellTask
from .ref import (CreateBWAIndices, CreateFASTAIndex, CreateSequenceDictionary,
                  FetchDbsnpVCF, FetchKnownIndelVCFs, FetchReferenceFASTA)
from .trim import TrimAdapters


@requires(TrimAdapters, FetchReferenceFASTA, CreateFASTAIndex,
          CreateBWAIndices)
class AlignReads(ShellTask):
    read_group = luigi.DictParameter()
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        '{0}.trim.{1}.{2}'.format(
                            parse_fq_id(fq_path=self.input()[0][0].path),
                            Path(self.input()[1].path).stem, s
                        )
                    )
                )
            ) for s in ['cram', 'cram.crai']
        ]

    def run(self):
        cram_path = self.output()[0].path
        run_id = Path(cram_path).stem
        print_log(f'Align reads:\t{run_id}')
        bwa = self.cf['bwa']
        samtools = self.cf['samtools']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_per_thread = self.cf['samtools_memory_per_thread']
        rg = '\\t'.join(
            [
                '@RG',
                'ID:{}'.format(self.read_group.get('ID') or 1),
                'PU:{}'.format(self.read_group.get('PU') or 'UNIT-1'),
                'SM:{}'.format(
                    self.read_group.get('SM')
                    or parse_fq_id(fq_path=self.input()[0][0].path)
                ),
                'PL:{}'.format(self.read_group.get('PL') or 'ILLUMINA'),
                'LB:{}'.format(self.read_group.get('LB') or 'LIBRARY-1')
            ] + [
                f'{k}:{v}' for k, v in self.read_group.items()
                if k not in ['ID', 'PU', 'SM', 'PL', 'LB']
            ]
        )
        fq_paths = [i.path for i in self.input()[0]]
        fa_path = self.input()[1].path
        index_paths = [o.path for o in self.input()[3]]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[bwa, samtools], cwd=self.cf['align_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                'set -eo pipefail && '
                + f'{bwa} mem -t {n_cpu} -R \'{rg}\' -Y {fa_path}'
                + ''.join([f' {a}' for a in fq_paths])
                + f' | {samtools} view -bS -'
                + f' | {samtools} sort -@ {n_cpu} -m {memory_per_thread}'
                + f' -T {cram_path}.sort -'
                + f' | {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                + f' -o {cram_path} -'
            ),
            input_files=[fa_path, *index_paths, *fq_paths],
            output_files=cram_path
        )
        self.run_shell(
            args=[
                f'set -e && {samtools} quickcheck -v {cram_path}',
                f'set -e && {samtools} index -@ {n_cpu} {cram_path}'
            ],
            input_files=cram_path, output_files=f'{cram_path}.crai'
        )


@requires(AlignReads, FetchReferenceFASTA, CreateFASTAIndex)
class MarkDuplicates(ShellTask):
    cf = luigi.DictParameter()
    priority = 70

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        Path(self.input()[0][0].path).stem + f'.markdup.{s}'
                    )
                )
            ) for s in ['cram', 'cram.crai', 'metrics.txt']
        ]

    def run(self):
        input_cram_path = self.input()[0][0].path
        run_id = Path(input_cram_path).stem
        print_log(f'Mark duplicates:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        samtools = self.cf['samtools']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_per_thread = self.cf['samtools_memory_per_thread']
        output_cram_path = self.output()[0].path
        markdup_metrics_txt_path = self.output()[2].path
        fa_path = self.input()[1].path
        tmp_bam_paths = [
            re.sub(r'\.cram$', s, output_cram_path)
            for s in ['.unfixed.unsorted.bam', '.unfixed.bam', '.bam']
        ]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, samtools], cwd=self.cf['align_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} MarkDuplicates'
                + f' --INPUT {input_cram_path}'
                + f' --REFERENCE_SEQUENCE {fa_path}'
                + f' --METRICS_FILE {markdup_metrics_txt_path}'
                + f' --OUTPUT {tmp_bam_paths[0]}'
                + ' --ASSUME_SORT_ORDER coordinate'
            ),
            input_files=[input_cram_path, fa_path],
            output_files=tmp_bam_paths[0]
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {samtools} sort -@ {n_cpu}'
                    + f' -m {memory_per_thread} -T {tmp_bam_paths[0]}.sort'
                    + f' -o {tmp_bam_paths[1]} {tmp_bam_paths[0]}'
                ),
                f'rm -f {tmp_bam_paths[0]}'
            ],
            input_files=tmp_bam_paths[0], output_files=tmp_bam_paths[1]
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {gatk}{gatk_opts} SetNmMdAndUqTags'
                    + f' --INPUT {tmp_bam_paths[1]}'
                    + f' --OUTPUT {tmp_bam_paths[2]}'
                    + f' --REFERENCE_SEQUENCE {fa_path}'
                ),
                f'rm -f {tmp_bam_paths[1]}'
            ],
            input_files=tmp_bam_paths[1], output_files=tmp_bam_paths[2]
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                    + f' -o {output_cram_path} {tmp_bam_paths[2]}'
                ),
                f'rm -f {tmp_bam_paths[2]}'
            ],
            input_files=[tmp_bam_paths[2], fa_path],
            output_files=output_cram_path
        )
        self.run_shell(
            args=[
                f'set -e && {samtools} quickcheck -v {output_cram_path}',
                f'set -e && {samtools} index -@ {n_cpu} {output_cram_path}'
            ],
            input_files=output_cram_path,
            output_files=f'{output_cram_path}.crai'
        )


@requires(MarkDuplicates, FetchReferenceFASTA, CreateFASTAIndex,
          CreateSequenceDictionary, FetchDbsnpVCF, FetchKnownIndelVCFs)
class ApplyBQSR(ShellTask):
    cf = luigi.DictParameter()
    priority = 80

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        Path(self.input()[0][0].path).stem + f'.bqsr.{s}'
                    )
                )
            ) for s in ['cram', 'cram.crai', 'data.csv']
        ]

    def run(self):
        input_cram_path = self.input()[0][0].path
        run_id = Path(input_cram_path).stem
        print_log(f'Apply Base Quality Score Recalibration:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        samtools = self.cf['samtools']
        n_cpu = self.cf['n_cpu_per_worker']
        output_cram_path = self.output()[0].path
        fa_path = self.input()[1].path
        fa_dict_path = self.input()[3].path
        known_site_vcf_gz_paths = [
            self.input()[4][0].path, *[o[0].path for o in self.input()[5]]
        ]
        bqsr_csv_path = self.output()[2].path
        tmp_bam_path = re.sub(r'\.cram$', '.bam', output_cram_path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, samtools], cwd=self.cf['align_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} BaseRecalibrator'
                + f' --input {input_cram_path}'
                + f' --reference {fa_path}'
                + f' --output {bqsr_csv_path}'
                + ' --use-original-qualities'
                + ''.join([
                    f' --known-sites {p}' for p in known_site_vcf_gz_paths
                ])
            ),
            input_files=[
                input_cram_path, fa_path, fa_dict_path,
                *known_site_vcf_gz_paths
            ],
            output_files=bqsr_csv_path
        )
        self.run_shell(
            args=(
                f'set -e && {gatk}{gatk_opts} ApplyBQSR'
                + f' --input {input_cram_path}'
                + f' --reference {fa_path}'
                + f' --bqsr-recal-file {bqsr_csv_path}'
                + f' --output {tmp_bam_path}'
                + ' --static-quantized-quals 10'
                + ' --static-quantized-quals 20'
                + ' --static-quantized-quals 30'
                + ' --add-output-sam-program-record'
                + ' --use-original-qualities'
                + ' --create-output-bam-index false'
            ),
            input_files=[input_cram_path, fa_path, bqsr_csv_path],
            output_files=tmp_bam_path
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                    + f' -o {output_cram_path} {tmp_bam_path}'
                ),
                f'rm -f {tmp_bam_path}'
            ],
            input_files=[tmp_bam_path, fa_path], output_files=output_cram_path
        )
        self.run_shell(
            args=[
                f'set -e && {samtools} quickcheck -v {output_cram_path}',
                f'set -e && {samtools} index -@ {n_cpu} {output_cram_path}'
            ],
            input_files=output_cram_path,
            output_files=f'{output_cram_path}.crai'
        )


@requires(ApplyBQSR, FetchReferenceFASTA)
class RemoveDuplicatesAndUnmapped(ShellTask):
    cf = luigi.DictParameter()
    priority = 90

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        Path(self.input()[0][0].path).stem + f'.prep.{s}'
                    )
                )
            ) for s in ['cram', 'cram.crai']
        ]

    def run(self):
        input_cram_path = self.input()[0][0].path
        run_id = Path(input_cram_path).stem
        print_log(f'Remove duplicates and unmapped reads:\t{run_id}')
        samtools = self.cf['samtools']
        n_cpu = self.cf['n_cpu_per_worker']
        output_cram_path = self.output()[0].path
        fa_path = self.input()[1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=samtools, cwd=self.cf['align_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -e && {samtools} view -@ {n_cpu} -T {fa_path}'
                + f' -F 1028 -CS -o {output_cram_path} {input_cram_path}'
            ),
            input_files=[input_cram_path, fa_path],
            output_files=output_cram_path
        )
        self.run_shell(
            args=[
                f'set -e && {samtools} quickcheck -v {output_cram_path}',
                f'set -e && {samtools} index -@ {n_cpu} {output_cram_path}'
            ],
            input_files=output_cram_path,
            output_files=f'{output_cram_path}.crai'
        )


class PrepareCRAMs(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    known_indel_vcf_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        return [
            RemoveDuplicatesAndUnmapped(
                fq_paths=f, read_group=r, sample_name=n,
                ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths, cf=self.cf
            ) for f, r, n
            in zip(self.fq_list, self.read_groups, self.sample_names)
        ]

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
