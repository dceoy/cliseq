#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import parse_fq_id, print_log
from .base import ShellTask
from .ref import (CreateBWAIndices, CreateFASTAIndex, CreateSequenceDictionary,
                  FetchKnownSiteVCFs, FetchReferenceFASTA)
from .trim import TrimAdapters


@requires(TrimAdapters, FetchReferenceFASTA, CreateFASTAIndex,
          CreateBWAIndices)
class AlignReads(ShellTask):
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
        r = '\'@RG\\tID:None\\tSM:None\\tPL:ILLUMINA\\tLB:None\''
        fq_paths = [i.path for i in self.input()[0]]
        fa_path = self.input()[1].path
        index_paths = [o.path for o in self.input()[3]]
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=self.cf['align_dir_path'], env={'REF_CACHE': '.ref_cache'}
        )
        self.run_bash(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'{samtools} 2>&1 | grep -e "Version:"',
                (
                    'set -eo pipefail && '
                    + f'{bwa} mem -t {n_cpu} -R {r} {fa_path} '
                    + ' '.join(fq_paths)
                    + f' | {samtools} view -bS -'
                    + f' | {samtools} sort -@ {n_cpu} -m {memory_per_thread}'
                    + f' -T {cram_path}.sort -'
                    + f' | {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                    + f' -o {cram_path} -'
                )
            ],
            input_files=[fa_path, *index_paths, *fq_paths],
            output_files=cram_path
        )
        self.run_bash(
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
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=self.cf['align_dir_path'], env={'REF_CACHE': '.ref_cache'}
        )
        self.run_bash(
            args=[
                f'{gatk} --version',
                f'{samtools} 2>&1 | grep -e "Version:"',
                (
                    f'set -e && {gatk}{gatk_opts} MarkDuplicates'
                    + f' --INPUT {input_cram_path}'
                    + f' --REFERENCE_SEQUENCE {fa_path}'
                    + f' --METRICS_FILE {markdup_metrics_txt_path}'
                    + f' --OUTPUT {tmp_bam_paths[0]}'
                    + ' --ASSUME_SORT_ORDER coordinate'
                ),
                (
                    f'set -e && {samtools} sort -@ {n_cpu}'
                    + f' -m {memory_per_thread} -T {tmp_bam_paths[0]}.sort'
                    + f' -o {tmp_bam_paths[1]} {tmp_bam_paths[0]}'
                ),
                f'rm -f {tmp_bam_paths[0]}',
                (
                    f'set -e && {gatk}{gatk_opts} SetNmMdAndUqTags'
                    + f' --INPUT {tmp_bam_paths[1]}'
                    + f' --OUTPUT {tmp_bam_paths[2]}'
                    + f' --REFERENCE_SEQUENCE {fa_path}'
                ),
                f'rm -f {tmp_bam_paths[1]}',
                (
                    f'set -e && {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                    + f' -o {output_cram_path} {tmp_bam_paths[2]}'
                ),
                f'rm -f {tmp_bam_paths[2]}'
            ],
            input_files=[input_cram_path, fa_path],
            output_files=output_cram_path
        )
        self.run_bash(
            args=[
                f'set -e && {samtools} quickcheck -v {output_cram_path}',
                f'set -e && {samtools} index -@ {n_cpu} {output_cram_path}'
            ],
            input_files=output_cram_path,
            output_files=f'{output_cram_path}.crai'
        )


@requires(MarkDuplicates, FetchReferenceFASTA, CreateFASTAIndex,
          CreateSequenceDictionary, FetchKnownSiteVCFs)
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
        known_site_vcf_gz_paths = [o[0].path for o in self.input()[4]]
        bqsr_csv_path = self.output()[2].path
        tmp_bam_path = re.sub(r'\.cram$', '.bam', output_cram_path)
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=self.cf['align_dir_path'], env={'REF_CACHE': '.ref_cache'}
        )
        self.run_bash(
            args=[
                f'{gatk} --version',
                f'{samtools} 2>&1 | grep -e "Version:"',
                (
                    f'set -e && {gatk}{gatk_opts} BaseRecalibrator'
                    + f' --input {input_cram_path}'
                    + f' --reference {fa_path}'
                    + f' --output {bqsr_csv_path}'
                    + ' --use-original-qualities'
                    + ''.join([
                        f' --known-sites {p}' for p in known_site_vcf_gz_paths
                    ])
                )
            ],
            input_files=[
                input_cram_path, fa_path, fa_dict_path,
                *known_site_vcf_gz_paths
            ],
            output_files=bqsr_csv_path
        )
        self.run_bash(
            args=[
                (
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
                (
                    f'set -e && {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                    + f' -o {output_cram_path} {tmp_bam_path}'
                ),
                f'rm -f {tmp_bam_path}'
            ],
            input_files=[input_cram_path, fa_path, bqsr_csv_path],
            output_files=output_cram_path
        )
        self.run_bash(
            args=[
                f'set -e && {samtools} quickcheck -v {output_cram_path}',
                f'set -e && {samtools} index -@ {n_cpu} {output_cram_path}'
            ],
            input_files=output_cram_path,
            output_files=f'{output_cram_path}.crai'
        )


@requires(ApplyBQSR)
class RemoveDuplicatesOrUnmapped(ShellTask):
    cf = luigi.DictParameter()
    priority = 90

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        Path(self.input()[0].path).stem + f'.filtered.{s}'
                    )
                )
            ) for s in ['cram', 'cram.crai']
        ]

    def run(self):
        input_cram_path = self.input()[0].path
        run_id = Path(input_cram_path).stem
        print_log(f'Remove duplicate :\t{run_id}')
        samtools = self.cf['samtools']
        n_cpu = self.cf['n_cpu_per_worker']
        output_cram_path = self.output()[0].path
        fa_path = self.input()[1].path
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=self.cf['align_dir_path'], env={'REF_CACHE': '.ref_cache'}
        )
        self.run_bash(
            args=(
                'set -e && {samtools} view -@ {n_cpu} -T {fa_path}'
                + f' -F 1028 -CS -o {output_cram_path} {input_cram_path}'
            ),
            input_files=[input_cram_path, fa_path],
            output_files=output_cram_path
        )
        self.run_bash(
            args=[
                f'set -e && {samtools} quickcheck -v {output_cram_path}',
                f'set -e && {samtools} index -@ {n_cpu} {output_cram_path}'
            ],
            input_files=output_cram_path,
            output_files=f'{output_cram_path}.crai'
        )


class PrepareCRAMs(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_dict = luigi.DictParameter()
    known_site_vcf_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        return [
            ApplyBQSR(
                fq_paths=v, ref_fa_paths=self.ref_fa_paths,
                known_site_vcf_paths=self.known_site_vcf_paths, cf=self.cf
            ) for v in self.fq_dict.values()
        ]


if __name__ == '__main__':
    luigi.run()
