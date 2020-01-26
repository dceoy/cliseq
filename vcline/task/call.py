#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id, print_log
from .align import PrepareCRAMs
from .base import ShellTask
from .ref import CreateFASTAIndex, FetchReferenceFASTA


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex)
class CallSomaticVariantsWithMutect2(ShellTask):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.cf['call_dir_path']).joinpath(
                    create_matched_id(*[i[0].path for i in self.input()[0]])
                    + '.mutect2.somatic.vcf'
                )
            )
        )

    def run(self):
        run_id = '.'.join(Path(self.output().path).stem.split('.')[:-2])
        print_log(f'Call somatic short variants with Mutect2:\t{run_id}')
        gatk = self.cf['gatk']
        gatk_opts = ' --java-options "{}"'.format(self.cf['gatk_java_options'])
        samtools = self.cf['samtools']
        save_memory = str(self.cf['memory_mb_per_worker'] < 16 * 1024).lower()
        cram_paths = [i[0].path for i in self.input()[0]]
        filtered_vcf_path = self.output().path
        unfiltered_vcf_path = re.sub(
            r'(\.vcf)$', '.unfiltered\\1', filtered_vcf_path
        )
        fa_path = self.input()[1].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[gatk, samtools], cwd=self.cf['call_dir_path'],
            env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=[
                (
                    f'set -e && {gatk}{gatk_opts} Mutect2'
                    + f' --reference {fa_path}'
                    + f' --input {cram_paths[0]}'
                    + f' --input {cram_paths[1]}'
                    + f' --output {unfiltered_vcf_path}'
                    + f' --disable-bam-index-caching {save_memory}'
                ),
                (
                    f'set -e && {gatk}{gatk_opts} FilterMutectCalls'
                    + f' --reference {fa_path}'
                    + f' --variant {unfiltered_vcf_path}'
                    + f' --output {filtered_vcf_path}'
                )
            ],
            input_files=[*cram_paths, fa_path],
            output_files=[filtered_vcf_path, unfiltered_vcf_path]
        )


class CallVariants(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    known_site_vcf_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 100

    def requires(self):
        return [
            CallSomaticVariantsWithMutect2(
                fq_list=self.fq_list, read_groups=self.read_groups,
                ref_fa_paths=self.ref_fa_paths,
                known_site_vcf_paths=self.known_site_vcf_paths, cf=self.cf
            )
        ]

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
