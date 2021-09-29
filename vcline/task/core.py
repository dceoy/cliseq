#!/usr/bin/env python

import re
from pathlib import Path

from ftarc.task.core import ShellTask


class VclineTask(ShellTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def generate_version_commands(commands):
        for c in ([commands] if isinstance(commands, str) else commands):
            n = Path(c).name
            if n in {'java', 'snpEff'} or n.endswith('.jar'):
                yield f'{c} -version'
            elif n in {'bwa', 'msisensor'}:
                yield f'{c} 2>&1 | grep -e "Program:" -e "Version:"'
            elif n == 'wget':
                yield f'{c} --version | head -1'
            elif n in {'bwa-mem2', 'cnvkit.py'}:
                yield f'{c} version'
            elif n == 'picard':
                yield f'{c} CreateSequenceDictionary --version'
            elif n == 'vep':
                yield f'{c} | grep -6 -e "Versions:"'
            else:
                yield f'{c} --version'

    @staticmethod
    def create_matched_id(tumor_name, normal_name):
        frags = [Path(n).stem.split('.') for n in [tumor_name, normal_name]]
        if frags[0][-1] != frags[1][-1]:
            somatic_id = '.'.join('.'.join(f) for f in frags)
        else:
            n_common = 0
            for i in range(1, min([len(f) for f in frags])):
                if frags[0][-i] == frags[1][-i]:
                    n_common += 1
                else:
                    break
            somatic_id = '.'.join(frags[0][:-n_common] + frags[1])
        return somatic_id

    @staticmethod
    def generate_gatk_java_options(n_cpu=1, memory_mb=4096):
        return ' '.join([
            '-Dsamjdk.compression_level=5',
            '-Dsamjdk.use_async_io_read_samtools=true',
            '-Dsamjdk.use_async_io_write_samtools=true',
            '-Dsamjdk.use_async_io_write_tribble=false',
            '-Xmx{}m'.format(int(memory_mb)), '-XX:+UseParallelGC',
            '-XX:ParallelGCThreads={}'.format(int(n_cpu))
        ])

    @classmethod
    def samtools_index(cls, sam_path, samtools='samtools', n_cpu=1):
        cls.run_shell(
            args=(
                f'set -e && {samtools} quickcheck -v {sam_path}'
                + f' && {samtools} index -@ {n_cpu} {sam_path}'
            ),
            input_files_or_dirs=sam_path,
            output_files_or_dirs=re.sub(
                r'\.(cr|b)am$', '.\\1am.\\1ai', str(sam_path)
            )
        )

    @classmethod
    def samtools_view(cls, input_sam_path, fa_path, output_sam_path,
                      samtools='samtools', n_cpu=1, add_args=None,
                      index_sam=False, remove_input=False):
        cls.run_shell(
            args=(
                f'set -e && {samtools} quickcheck -v {input_sam_path}'
                + f' && {samtools} view -@ {n_cpu} -T {fa_path}'
                + ' -{0}S{1}'.format(
                    ('C' if str(output_sam_path).endswith('.cram') else 'b'),
                    (f' {add_args}' if add_args else '')
                )
                + f' -o {output_sam_path} {input_sam_path}'
            ),
            input_files_or_dirs=[
                input_sam_path, fa_path, f'{fa_path}.fai'
            ],
            output_files_or_dirs=output_sam_path
        )
        if index_sam:
            cls.samtools_index(
                sam_path=output_sam_path, samtools=samtools, n_cpu=n_cpu
            )
        if remove_input:
            cls.remove_files_and_dirs(input_sam_path)

    @classmethod
    def samtools_merge(cls, input_sam_paths, fa_path, output_sam_path,
                       samtools='samtools', n_cpu=1, memory_mb=1024,
                       index_sam=True, remove_input=True):
        memory_mb_per_thread = int(memory_mb / n_cpu / 8)
        cls.run_shell(
            args=(
                f'set -eo pipefail && {samtools} merge -@ {n_cpu} -r -'
                + ''.join(f' {p}' for p in input_sam_paths)
                + f' | {samtools} sort -@ {n_cpu} -m {memory_mb_per_thread}M'
                + f' -O bam -l 0 -T {output_sam_path}.sort -'
                + f' | {samtools} view -@ {n_cpu}'
                + f' -T {fa_path}'
                + ' -{}S'.format(
                    'C' if str(output_sam_path).endswith('.cram') else 'b'
                )
                + f' -o {output_sam_path} -'
            ),
            input_files_or_dirs=[
                *input_sam_paths, fa_path, f'{fa_path}.fai'
            ],
            output_files_or_dirs=output_sam_path
        )
        if index_sam:
            cls.samtools_index(
                sam_path=output_sam_path, samtools=samtools, n_cpu=n_cpu
            )
        if remove_input:
            cls.remove_files_and_dirs(*input_sam_paths)

    @classmethod
    def tabix_tbi(cls, tsv_path, tabix='tabix', preset='vcf', **kwargs):
        cls.run_shell(
            args=f'set -e && {tabix} --preset {preset} {tsv_path}',
            input_files_or_dirs=tsv_path,
            output_files_or_dirs=f'{tsv_path}.tbi', **kwargs
        )

    @classmethod
    def bcftools_index(cls, vcf_path, bcftools='bcftools', n_cpu=1, tbi=True):
        cls.run_shell(
            args=(
                f'set -e && {bcftools} index --threads {n_cpu}'
                + (' --tbi' if tbi else ' --csi')
                + f' {vcf_path}'
            ),
            input_files_or_dirs=vcf_path,
            output_files_or_dirs=re.sub(
                r'\.(bcf|vcf.gz)$', '.\\1.{}'.format('tbi' if tbi else 'csi'),
                str(vcf_path)
            )
        )

    @classmethod
    def bcftools_concat(cls, input_vcf_paths, output_vcf_path,
                        bcftools='bcftools', n_cpu=1, memory_mb=1024,
                        index_vcf=True, remove_input=True):
        cls.run_shell(
            args=(
                f'set -eo pipefail && {bcftools} concat --threads {n_cpu}'
                + ''.join(f' {p}' for p in input_vcf_paths)
                + f' | {bcftools} sort --max-mem {memory_mb}M'
                + f' --temp-dir {output_vcf_path}.sort --output-type z'
                + f' --output-file {output_vcf_path} -'
            ),
            input_files_or_dirs=input_vcf_paths,
            output_files_or_dirs=output_vcf_path
        )
        if index_vcf:
            cls.bcftools_index(
                vcf_path=output_vcf_path, bcftools=bcftools, n_cpu=n_cpu,
                tbi=True
            )
        if remove_input:
            cls.remove_files_and_dirs(*input_vcf_paths)

    @classmethod
    def bcftools_sort(cls, input_vcf_path, output_vcf_path,
                      bcftools='bcftools', n_cpu=1, memory_mb=1024,
                      index_vcf=True, remove_input=True):
        cls.run_shell(
            args=(
                f'set -e && {bcftools} sort --max-mem {memory_mb}M'
                + f' --temp-dir {output_vcf_path}.sort'
                + f' --output-type z --output-file {output_vcf_path}'
                + f' {input_vcf_path}'
            ),
            input_files_or_dirs=input_vcf_path,
            output_files_or_dirs=output_vcf_path
        )
        if index_vcf:
            cls.bcftools_index(
                vcf_path=output_vcf_path, bcftools=bcftools, n_cpu=n_cpu,
                tbi=True
            )
        if remove_input:
            cls.remove_files_and_dirs(input_vcf_path)

    @classmethod
    def picard_mergevcfs(cls, input_vcf_paths, output_vcf_path,
                         picard='picard', remove_input=True):
        cls.run_shell(
            args=(
                f'set -e && {picard} MergeVcfs'
                + ''.join(f' --INPUT {v}' for v in input_vcf_paths)
                + f' --OUTPUT {output_vcf_path}'
            ),
            input_files_or_dirs=input_vcf_paths,
            output_files_or_dirs=[output_vcf_path, f'{output_vcf_path}.tbi']
        )
        if remove_input:
            cls.remove_files_and_dirs(*input_vcf_paths)
