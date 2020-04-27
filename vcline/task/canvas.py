#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareBAMNormal, PrepareBAMTumor
from .base import ShellTask
from .haplotypecaller import GenotypeHaplotypeCallerGVCFVCF
from .mutect2 import CallVariantsWithMutect2
from .ref import (FetchReferenceFASTA, FetchResourceFASTA, FetchResourceFile,
                  UncompressCnvBlackListBED)


class CreateCanvasGenomeSymlinks(ShellTask):
    ref_fa_path = luigi.Parameter()
    genomesize_xml_path = luigi.Parameter()
    kmer_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def requires(self):
        return [
            FetchReferenceFASTA(ref_fa_path=self.ref_fa_path, cf=self.cf),
            FetchResourceFile(
                resource_file_path=self.genomesize_xml_path, cf=self.cf
            ),
            FetchResourceFASTA(
                resource_file_path=self.kmer_fa_path, cf=self.cf
            )
        ]

    def output(self):
        return [
            luigi.LocalTarget(
                str(Path(self.cf['ref_dir_path']).joinpath(f'canvas/{n}'))
            ) for n in [
                'genome.fa', 'genome.fa.fai', 'GenomeSize.xml',
                'kmer.fa', 'kmer.fa.fai'
            ]
        ]

    def run(self):
        src_paths = [
            *[i.path for i in self.input()[0]], self.input()[1].path,
            *[i.path for i in self.input()[2]]
        ]
        run_id = Path(src_paths[0]).stem
        self.print_log(f'Prepare Canvas genome folder:\t{run_id}')
        symlinks = {
            '../{}'.format(Path(p).name): o.path
            for p, o in zip(src_paths, self.output())
        }
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=str(Path(list(symlinks.values())[0]).parent),
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=[f'ln -s {s} {d}' for s, d in symlinks.items()],
            output_files_or_dirs=list(symlinks.values())
        )


@requires(PrepareBAMTumor, PrepareBAMNormal, CreateCanvasGenomeSymlinks,
          UncompressCnvBlackListBED, GenotypeHaplotypeCallerGVCFVCF,
          CallVariantsWithMutect2)
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
        bam_path = self.input()[0][0].path
        run_id = Path(output_dir_path).name
        self.print_log(f'Call somatic CNVs with Canvas:\t{run_id}')
        canvas = self.cf['Canvas']
        sample_name = self.sample_names[0]
        kmer_fa_path = self.input()[2][3].path
        canvas_genome_dir_path = str(Path(self.input()[2][0].path).parent)
        filter_bed_path = self.input()[3].path
        germline_b_allele_vcf = self.input()[4][0].path
        somatic_vcf = self.input()[5][0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=canvas, cwd=self.cf['somatic_cnv_canvas_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {canvas} Somatic-WGS'
                + f' --bam={bam_path}'
                + f' --sample-name {sample_name}'
                + f' --reference {kmer_fa_path}'
                + f' --genome-folder {canvas_genome_dir_path}'
                + f' --filter-bed {filter_bed_path}'
                + f' --output {output_dir_path}'
                + f' --sample-b-allele-vcf {germline_b_allele_vcf}'
                + f' --somatic-vcf {somatic_vcf}'
            ),
            input_files_or_dirs=[
                bam_path, kmer_fa_path, canvas_genome_dir_path,
                filter_bed_path
            ],
            output_files_or_dirs=output_dir_path
        )


if __name__ == '__main__':
    luigi.run()
