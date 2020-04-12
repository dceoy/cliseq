#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareTumorCRAM
from .base import ShellTask
from .ref import (CreateCnvBlackListBED, FlagUniqueKmers,
                  PrepareCanvasGenomeFolder)
from .strelka import (CallGermlineVariantsWithStrelka,
                      CallSomaticVariantsWithStrelka)


@requires(PrepareTumorCRAM, FlagUniqueKmers, PrepareCanvasGenomeFolder,
          CreateCnvBlackListBED, CallGermlineVariantsWithStrelka,
          CallSomaticVariantsWithStrelka)
class CallSomaticCopyNumberVariantsWithCanvas(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['somatic_cnv_canvas_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0:2]]
                        )
                    )
                )
            )
        ]

    def run(self):
        output_dir_path = self.output()[0].path
        run_id = Path(output_dir_path).name
        self.print_log(f'Call somatic CNVs with Canvas:\t{run_id}')
        canvas = self.cf['Canvas']
        sample_name = self.sample_names[0]
        input_cram_path = self.input()[0][0].path
        fa_path = self.input()[1][0].path
        fai_path = self.input()[1][1].path
        canvas_genome_dir_path = self.input()[2].path
        filter_bed_path = self.input()[3][0].path
        germline_b_allele_vcf = [
            i.path for i in self.input()[4]
            if i.path.endswith('.variants.vcf.gz')
        ][0]
        somatic_vcf = [
            i.path for i in self.input()[5]
            if i.path.endswith('.somatic.vcf.gz')
        ][0]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=canvas, cwd=self.cf['somatic_cnv_canvas_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {canvas} Somatic-WGS'
                + f' --bam={input_cram_path}'
                + f' --sample-name {sample_name}'
                + f' --reference {fa_path}'
                + f' --genome-folder {canvas_genome_dir_path}'
                + f' --filter-bed {filter_bed_path}'
                + f' --output {output_dir_path}'
                + f' --sample-b-allele-vcf={germline_b_allele_vcf}'
                + f' --somatic-vcf {somatic_vcf}'
            ),
            input_files_or_dirs=[
                input_cram_path, fa_path, fai_path, canvas_genome_dir_path,
                filter_bed_path
            ],
            output_files_or_dirs=output_dir_path
        )


if __name__ == '__main__':
    luigi.run()
