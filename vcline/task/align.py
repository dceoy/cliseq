#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import parse_fq_id, print_log
from .base import ShellTask
from .ref import CreateBWAIndexes, CreateFASTAIndex, FetchGenomeFASTA
from .trim import TrimAdapters


@requires(FetchGenomeFASTA, CreateBWAIndexes, CreateFASTAIndex, TrimAdapters)
class AlignReads(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        '{0}.trim.{1}.cram{2}'.format(
                            parse_fq_id(fq_path=self.input()[3][0].path),
                            Path(self.input()[0].path).stem, s
                        )
                    )
                )
            ) for s in ['', '.crai']
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
        fa_path = self.input()[0].path
        index_paths = [o.path for o in [*self.input()[1], self.input()[2]]]
        fq_paths = [i.path for i in self.input()[3]]
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['align_dir_path']
        )
        self.run_bash(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'{samtools} 2>&1 | grep -e "Version:"',
                (
                    'set -eo pipefail && '
                    + f'{bwa} mem -t {n_cpu} -R {r} {fa_path} '
                    + ' '.join(fq_paths)
                    + f' | {samtools} view -@ {n_cpu} -bS -o - -'
                    + f' | {samtools} sort -@ {n_cpu} -m {memory_per_thread}'
                    + f' -T {cram_path}.sort -o - -'
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


@requires(FetchGenomeFASTA, CreateFASTAIndex, AlignReads)
class MarkDuplicates(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['align_dir_path']).joinpath(
                        '{0}.markdup.cram{1}'.format(
                            str(Path(self.input()[2][0].path).stem), s
                        )
                    )
                )
            ) for s in ['', '.crai']
        ]

    def run(self):
        input_cram_path = self.input()[2][0].path
        run_id = Path(input_cram_path).stem
        print_log(f'Mark duplicates:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        samtools = self.cf['samtools']
        n_cpu = self.cf['n_cpu_per_worker']
        output_cram_path = self.output()[0].path
        fa_path = self.input()[0].path
        fai_path = self.input()[1].path
        markdup_metrics_txt_path = str(
            Path(output_cram_path).parent.joinpath(
                Path(output_cram_path).stem + '.metrics.txt'
            )
        )
        tmp_bam_path = f'{output_cram_path}.tmp.bam'
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['align_dir_path']
        )
        self.run_bash(
            args=[
                f'{gatk} --version',
                f'{samtools} 2>&1 | grep -e "Version:"',
                (
                    'set -eo pipefail && '
                    + f'{gatk}{gatk_opts} MarkDuplicates'
                    + f' --INPUT {input_cram_path}'
                    + f' --REFERENCE_SEQUENCE {fa_path}'
                    + f' --METRICS_FILE {markdup_metrics_txt_path}'
                    + f' --OUTPUT {tmp_bam_path}'
                    + ' --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500'
                    + ' --ASSUME_SORT_ORDER coordinate'
                )
            ],
            input_files=[input_cram_path, fa_path, fai_path],
            output_files=tmp_bam_path
        )
        self.run_bash(
            args=[
                f'set -e && {samtools} quickcheck -v {tmp_bam_path}',
                (
                    'set -eo pipefail && '
                    + f'{gatk}{gatk_opts} SortSam'
                    + f'--INPUT {tmp_bam_path}'
                    + ' --OUTPUT /dev/stdout'
                    + ' --SORT_ORDER coordinate'
                    + f' | {gatk}{gatk_opts} SetNmMdAndUqTags'
                    + ' --INPUT /dev/stdin'
                    + ' --OUTPUT /dev/stdout'
                    + f' --REFERENCE_SEQUENCE {fa_path}'
                    + f' | {samtools} view -@ {n_cpu} -T {fa_path} -CS'
                    + f' -o {output_cram_path} -'
                ),
                f'set -e && rm -f {tmp_bam_path}'
            ],
            input_files=[tmp_bam_path, fa_path, fai_path],
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


@requires(FetchGenomeFASTA, CreateFASTAIndex, MarkDuplicates)
class ApplyBQSR(luigi.WrapperTask):
    priority = 10


if __name__ == '__main__':
    luigi.run()
