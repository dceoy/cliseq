#!/usr/bin/env python

import sys
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareBAMNormal, PrepareBAMTumor
from .base import ShellTask
from .haplotypecaller import GenotypeHaplotypeCallerGVCF
from .mutect2 import CallVariantsWithMutect2
from .ref import (CreateCnvBlackListBED, FetchReferenceFASTA,
                  FetchResourceFASTA, FetchResourceFile, UncompressBgzipFiles)


class CreateCanvasGenomeSymlinks(ShellTask):
    ref_fa_path = luigi.Parameter()
    genomesize_xml_path = luigi.Parameter()
    kmer_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 50

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
                Path(self.cf['ref_dir_path']).joinpath(f'canvas/{n}')
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
            cwd=Path(list(symlinks.values())[0]).parent,
            remove_if_failed=self.cf['remove_if_failed']
        )
        for s, d in symlinks.items():
            self.run_shell(args=f'ln -s {s} {d}', output_files_or_dirs=d)


@requires(CreateCnvBlackListBED)
class UncompressCnvBlackListBED(luigi.Task):
    cf = luigi.DictParameter()
    priority = 50

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath(
                'canvas/{}'.format(Path(self.input()[0].path).stem)
            )
        )

    def run(self):
        yield UncompressBgzipFiles(
            bgz_paths=[self.input()[0].path],
            dest_dir_path=Path(self.output().path).parent, cf=self.cf
        )


class CreateUniqueRegionManifest(ShellTask):
    exome_manifest_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 10

    def requires(self):
        return FetchResourceFile(
            resource_file_path=self.exome_manifest_path, cf=self.cf
        )

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['somatic_cnv_canvas_dir_path']).joinpath(
                Path(self.input().path).stem + '.uniq.txt'
            )
        )

    def run(self):
        output_manifest_path = self.output().path
        run_id = Path(output_manifest_path).name
        self.print_log(f'Create a unique region maifest:\t{run_id}')
        input_manifest_path = self.input().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=self.cf['somatic_cnv_canvas_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && {sys.executable}'
                + ' -c \'{}\''.format(
                    'from fileinput import input; '
                    'import sys; '
                    '[sys.stdout.write('
                    '"{0}_{1}\t{1}\t{2}".format(*s.split("\t", maxsplit=2))'
                    ' if i > 1 else s)'
                    ' for i, s in enumerate(input())];'
                ) + f' {input_manifest_path} > {output_manifest_path}'
            ),
            input_files_or_dirs=input_manifest_path,
            output_files_or_dirs=output_manifest_path
        )


@requires(PrepareBAMTumor, PrepareBAMNormal, CreateCanvasGenomeSymlinks,
          UncompressCnvBlackListBED, GenotypeHaplotypeCallerGVCF,
          CallVariantsWithMutect2)
class CallSomaticCopyNumberVariantsWithCanvas(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    exome_manifest_path = luigi.Parameter(default='')
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['somatic_cnv_canvas_dir_path']).joinpath(
                    create_matched_id(*[i[0].path for i in self.input()[0:2]])
                )
            )
        ]

    def run(self):
        output_dir_path = self.output()[0].path
        tumor_bam_path = self.input()[0][0].path
        run_id = Path(output_dir_path).name
        self.print_log(f'Call somatic CNVs with Canvas:\t{run_id}')
        canvas = self.cf['Canvas']
        sample_name = self.sample_names[0]
        kmer_fa_path = self.input()[2][3].path
        canvas_genome_dir_path = str(Path(self.input()[2][0].path).parent)
        filter_bed_path = self.input()[3].path
        germline_b_allele_vcf_path = self.input()[4][0].path
        somatic_vcf_path = self.input()[5][0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=canvas, cwd=self.cf['somatic_cnv_canvas_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        if self.cf['exome']:
            target = yield CreateUniqueRegionManifest(
                exome_manifest_path=self.exome_manifest_path, cf=self.cf
            )
            manifest_path = target.path
            self.run_shell(
                args=(
                    f'set -e && {canvas} Somatic-Enrichment'
                    + f' --bam={tumor_bam_path}'
                    + f' --manifest={manifest_path}'
                    + f' --sample-b-allele-vcf {germline_b_allele_vcf_path}'
                    + f' --sample-name {sample_name}'
                    + f' --reference {kmer_fa_path}'
                    + f' --genome-folder {canvas_genome_dir_path}'
                    + f' --filter-bed {filter_bed_path}'
                    + f' --output {output_dir_path}'
                ),
                input_files_or_dirs=[
                    tumor_bam_path, somatic_vcf_path,
                    germline_b_allele_vcf_path, manifest_path, kmer_fa_path,
                    canvas_genome_dir_path, filter_bed_path
                ],
                output_files_or_dirs=output_dir_path
            )
        else:
            self.run_shell(
                args=(
                    f'set -e && {canvas} Somatic-WGS'
                    + f' --bam={tumor_bam_path}'
                    + f' --somatic-vcf {somatic_vcf_path}'
                    + f' --sample-b-allele-vcf {germline_b_allele_vcf_path}'
                    + f' --sample-name {sample_name}'
                    + f' --reference {kmer_fa_path}'
                    + f' --genome-folder {canvas_genome_dir_path}'
                    + f' --filter-bed {filter_bed_path}'
                    + f' --output {output_dir_path}'
                ),
                input_files_or_dirs=[
                    tumor_bam_path, somatic_vcf_path,
                    germline_b_allele_vcf_path, kmer_fa_path,
                    canvas_genome_dir_path, filter_bed_path
                ],
                output_files_or_dirs=output_dir_path
            )


if __name__ == '__main__':
    luigi.run()
