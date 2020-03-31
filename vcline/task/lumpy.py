#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMs
from .base import ShellTask
from .ref import CreateFASTAIndex, FetchReferenceFASTA


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex)
class CallStructualVariantsWithLumpy(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['lumpy_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0]]
                        ) + f'.lumpy.{s}'
                    )
                )
            ) for s in ['vcf.gz', 'vcf.gz.tbi']
        ]

    def run(self):
        output_vcf_path = self.output()[0].path
        run_id = Path(Path(output_vcf_path).stem).stem
        self.print_log(f'Call somatic SVs with Lumpy:\t{run_id}')
        lumpy = self.cf['lumpy']
        lumpyexpress = self.cf['lumpyexpress']
        extractsplitreads_bwamem = str(
            Path(self.cf['lumpyexpress']).parent.parent.joinpath(
                'scripts/extractSplitReads_BwaMem'
            )
        )
        python = self.cf['python']
        samblaster = self.cf['samblaster']
        sambamba = self.cf['sambamba']
        samtools = self.cf['samtools']
        gawk = self.cf['gawk']
        bgzip = self.cf['bgzip']
        bcftools = self.cf['bcftools']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_per_thread = self.cf['samtools_memory_per_thread']
        input_cram_paths = [i[0].path for i in self.input()[0]]
        fa_path = self.input()[1].path
        fai_path = self.input()[2].path
        discordants_bam_paths = [
            str(
                Path(self.cf['lumpy_dir_path']).joinpath(
                    Path(c).stem + '.discordants.bam'
                )
            ) for c in input_cram_paths
        ]
        splitters_bam_paths = [
            str(
                Path(self.cf['lumpy_dir_path']).joinpath(
                    Path(c).stem + '.splitters.bam'
                )
            ) for c in input_cram_paths
        ]
        uncompressed_vcf_path = re.sub(r'\.gz', '', output_vcf_path)
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[
                lumpy, python, samblaster, sambamba, samtools, gawk, bgzip,
                bcftools
            ],
            cwd=self.cf['lumpy_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=[
                (
                    f'set -eo pipefail && '
                    + f'{samtools} view -@ {n_cpu} -T {fa_path} -bS'
                    + f' -F 1294 -o {o} {i}'
                ) for i, o in zip(input_cram_paths, discordants_bam_paths)
            ],
            input_files_or_dirs=[*input_cram_paths, fa_path, fai_path],
            output_files_or_dirs=discordants_bam_paths
        )
        self.run_shell(
            args=[
                (
                    f'set -eo pipefail && '
                    + f'{samtools} view -@ {n_cpu} -T {fa_path} -S -h {i}'
                    + f' | {extractsplitreads_bwamem} -i stdin'
                    + f' | {samtools} view -@ {n_cpu} -bS -'
                    + f' | {samtools} sort -@ {n_cpu} -m {memory_per_thread}'
                    + f' -T {o}.sort -o {o} -'
                ) for i, o in zip(input_cram_paths, splitters_bam_paths)
            ],
            input_files_or_dirs=[*input_cram_paths, fa_path, fai_path],
            output_files_or_dirs=splitters_bam_paths
        )
        self.run_shell(
            args=(
                f'set -e && {lumpyexpress}'
                + ' -B {0},{1}'.format(*input_cram_paths)
                + ' -S {0},{1}'.format(*splitters_bam_paths)
                + ' -D {0},{1}'.format(*discordants_bam_paths)
                + f' -R {fa_path}'
                + f' -T {uncompressed_vcf_path}.temp'
                + f' -t {n_cpu}'
                + f' -o {uncompressed_vcf_path}'
            ),
            input_files_or_dirs=[*input_cram_paths, fa_path, fai_path],
            output_files_or_dirs=uncompressed_vcf_path
        )
        self.run_shell(
            args=f'set -e && {bgzip} -@ {n_cpu} {uncompressed_vcf_path}',
            input_files_or_dirs=uncompressed_vcf_path,
            output_files_or_dirs=output_vcf_path
        )
        self.run_shell(
            args=(
                f'set -e && {bcftools} index'
                + ' --tbi'
                + f' --threads {n_cpu}'
                + f' {output_vcf_path}'
            ),
            input_files_or_dirs=output_vcf_path,
            output_files_or_dirs=f'{output_vcf_path}.tbi'
        )


if __name__ == '__main__':
    luigi.run()
