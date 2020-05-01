#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import parse_fq_id
from .base import BaseTask, ShellTask
from .ref import (CreateBWAIndices, CreateSequenceDictionary, FetchDbsnpVCF,
                  FetchKnownIndelVCF, FetchMillsIndelVCF, FetchReferenceFASTA)
from .samtools import SamtoolsIndex, SamtoolsView, SortSAM
from .trim import TrimAdapters


@requires(TrimAdapters, FetchReferenceFASTA, CreateBWAIndices)
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
                            Path(self.input()[1][0].path).stem, s
                        )
                    )
                )
            ) for s in ['cram', 'cram.crai']
        ]

    def run(self):
        cram_path = self.output()[0].path
        run_id = Path(cram_path).stem
        self.print_log(f'Align reads:\t{run_id}')
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
        fa_path = self.input()[1][0].path
        index_paths = [o.path for o in self.input()[2]]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[bwa, samtools], cwd=self.cf['align_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                'set -eo pipefail && '
                + f'{bwa} mem -t {n_cpu} -R \'{rg}\' -T 0 -P {fa_path}'
                + ''.join([f' {a}' for a in fq_paths])
                + f' | {samtools} sort -@ {n_cpu} -m {memory_per_thread}'
                + f' -O bam -l 0 -T {cram_path}.sort -'
                + f' | {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                + f' -o {cram_path} -'
            ),
            input_files_or_dirs=[fa_path, *index_paths, *fq_paths],
            output_files_or_dirs=cram_path
        )
        yield SamtoolsIndex(
            sam_path=cram_path, samtools=samtools, n_cpu=n_cpu,
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )


@requires(AlignReads, FetchReferenceFASTA, CreateSequenceDictionary)
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
        self.print_log(f'Mark duplicates:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        output_cram_path = self.output()[0].path
        markdup_metrics_txt_path = self.output()[2].path
        fa_path = self.input()[1][0].path
        tmp_bam_paths = [
            re.sub(r'\.cram$', f'.{s}.bam', output_cram_path)
            for s in ['unfixed.unsorted', 'unfixed']
        ]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['align_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
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
            input_files_or_dirs=[input_cram_path, fa_path],
            output_files_or_dirs=tmp_bam_paths[0]
        )
        yield SortSAM(
            input_sam_path=tmp_bam_paths[0],
            output_sam_path=tmp_bam_paths[1], fa_path=fa_path,
            samtools=self.cf['samtools'], n_cpu=self.cf['n_cpu_per_worker'],
            memory_per_thread=self.cf['samtools_memory_per_thread'],
            remove_input=True, index_sam=False,
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        yield SetNmMdAndUqTags(
            input_bam_path=tmp_bam_paths[1], fa_path=fa_path,
            output_cram_path=output_cram_path, gatk=gatk, gatk_opts=gatk_opts,
            samtools=self.cf['samtools'], n_cpu=self.cf['n_cpu_per_worker'],
            memory_per_thread=self.cf['samtools_memory_per_thread'],
            remove_input=True, log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )


class SetNmMdAndUqTags(ShellTask):
    input_bam_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    output_cram_path = luigi.Parameter()
    gatk = luigi.Parameter()
    gatk_opts = luigi.Parameter()
    samtools = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_per_thread = luigi.Parameter(default='768M')
    remove_input = luigi.BoolParameter(default=True)
    log_dir_path = luigi.Parameter(default='')
    remove_if_failed = luigi.BoolParameter(default=True)
    priority = 70

    def output(self):
        return [
            luigi.LocalTarget(self.output_cram_path + s) for s in ['', '.crai']
        ]

    def run(self):
        run_id = Path(self.input_bam_path).stem
        self.print_log(f'Calculatesthe NM, MD, and UQ tags:\t{run_id}')
        tmp_bam_path = re.sub(r'\.cram$', '.bam', self.output_cram_path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.log_dir_path,
            commands=self.gatk, cwd=str(Path(tmp_bam_path).parent),
            remove_if_failed=self.remove_if_failed,
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk}{self.gatk_opts} SetNmMdAndUqTags'
                + f' --INPUT {self.input_bam_path}'
                + f' --OUTPUT {tmp_bam_path}'
                + f' --REFERENCE_SEQUENCE {self.fa_path}'
            ),
            input_files_or_dirs=[self.input_bam_path, self.fa_path],
            output_files_or_dirs=tmp_bam_path
        )
        yield SamtoolsView(
            input_sam_path=tmp_bam_path,
            output_sam_path=self.output_cram_path,
            fa_path=self.fa_path, samtools=self.samtools, n_cpu=self.n_cpu,
            remove_input=self.remove_input, log_dir_path=self.log_dir_path,
            remove_if_failed=self.remove_if_failed
        )
        yield SamtoolsIndex(
            sam_path=self.output_cram_path, samtools=self.samtools,
            n_cpu=self.n_cpu, log_dir_path=self.log_dir_path,
            remove_if_failed=self.remove_if_failed
        )


@requires(MarkDuplicates, FetchReferenceFASTA, CreateSequenceDictionary,
          FetchDbsnpVCF, FetchMillsIndelVCF, FetchKnownIndelVCF)
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
        self.print_log(f'Apply Base Quality Score Recalibration:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        output_cram_path = self.output()[0].path
        fa_path = self.input()[1][0].path
        fa_dict_path = self.input()[2].path
        known_site_vcf_gz_paths = [i[0].path for i in self.input()[3:6]]
        bqsr_csv_path = self.output()[2].path
        tmp_bam_path = re.sub(r'\.cram$', '.bam', output_cram_path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=gatk, cwd=self.cf['align_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
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
            input_files_or_dirs=[
                input_cram_path, fa_path, fa_dict_path,
                *known_site_vcf_gz_paths
            ],
            output_files_or_dirs=bqsr_csv_path
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
            input_files_or_dirs=[input_cram_path, fa_path, bqsr_csv_path],
            output_files_or_dirs=tmp_bam_path
        )
        yield SamtoolsView(
            input_sam_path=tmp_bam_path, output_sam_path=output_cram_path,
            fa_path=fa_path, samtools=self.cf['samtools'],
            n_cpu=self.cf['n_cpu_per_worker'], remove_input=True,
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        yield SamtoolsIndex(
            sam_path=output_cram_path, samtools=self.cf['samtools'],
            n_cpu=self.cf['n_cpu_per_worker'],
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )


@requires(ApplyBQSR, FetchReferenceFASTA)
class RemoveDuplicates(BaseTask):
    cf = luigi.DictParameter()
    priority = 90

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        Path(self.input()[0][0].path).stem + f'.dedup.{s}'
                    )
                )
            ) for s in ['cram', 'cram.crai']
        ]

    def run(self):
        yield SamtoolsView(
            input_sam_path=self.input()[0][0].path,
            output_sam_path=self.output()[0].path,
            fa_path=self.input()[1][0].path, samtools=self.cf['samtools'],
            n_cpu=self.cf['n_cpu_per_worker'], add_args='-F 1024',
            message='Remove duplicates', remove_input=False,
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        yield SamtoolsIndex(
            sam_path=self.output()[0].path, samtools=self.cf['samtools'],
            n_cpu=self.cf['n_cpu_per_worker'],
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )


class PrepareCRAMTumor(luigi.WrapperTask):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        if self.cram_list:
            return super().requires()
        else:
            return RemoveDuplicates(
                fq_paths=self.fq_list[0], read_group=self.read_groups[0],
                sample_name=self.sample_names[0], ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path, cf=self.cf
            )

    def output(self):
        if self.cram_list:
            return [
                luigi.LocalTarget(self.cram_list[0] + s) for s in ['', '.crai']
            ]
        else:
            return self.input()


class PrepareCRAMNormal(luigi.WrapperTask):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        if self.cram_list:
            return super().requires()
        else:
            return RemoveDuplicates(
                fq_paths=self.fq_list[1], read_group=self.read_groups[1],
                sample_name=self.sample_names[1], ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path, cf=self.cf
            )

    def output(self):
        if self.cram_list:
            return [
                luigi.LocalTarget(self.cram_list[1] + s) for s in ['', '.crai']
            ]
        else:
            return self.input()


@requires(PrepareCRAMTumor, PrepareCRAMNormal)
class PrepareCRAMsMatched(luigi.WrapperTask):
    priority = 100

    def output(self):
        return self.input()


@requires(PrepareCRAMTumor, FetchReferenceFASTA)
class PrepareBAMTumor(BaseTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        Path(self.input()[0][0].path).stem + f'.bam{s}'
                    )
                )
            ) for s in ['', '.bai']
        ]

    def run(self):
        cram_path = self.input()[0][0].path
        bam_path = self.output()[0].path
        yield SamtoolsView(
            input_sam_path=cram_path, output_sam_path=bam_path,
            fa_path=self.input()[1][0].path, samtools=self.cf['samtools'],
            n_cpu=self.cf['n_cpu_per_worker'], message='Convert CRAM to BAM',
            remove_input=False, log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
        )
        yield SamtoolsIndex(
            sam_path=bam_path, samtools=self.cf['samtools'],
            n_cpu=self.cf['n_cpu_per_worker'],
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )


@requires(PrepareCRAMNormal, FetchReferenceFASTA)
class PrepareBAMNormal(BaseTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        Path(self.input()[0][0].path).stem + f'.bam{s}'
                    )
                )
            ) for s in ['', 'bai']
        ]

    def run(self):
        cram_path = self.input()[0][0].path
        bam_path = self.output()[0].path
        yield SamtoolsView(
            input_sam_path=cram_path, output_sam_path=bam_path,
            fa_path=self.input()[1][0].path, samtools=self.cf['samtools'],
            n_cpu=self.cf['n_cpu_per_worker'], message='Convert CRAM to BAM',
            remove_input=False, log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
        )
        yield SamtoolsIndex(
            sam_path=bam_path, samtools=self.cf['samtools'],
            n_cpu=self.cf['n_cpu_per_worker'],
            log_dir_path=self.cf['log_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )


if __name__ == '__main__':
    luigi.run()
