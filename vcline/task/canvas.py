#!/usr/bin/env python

import sys
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .base import ShellTask
from .haplotypecaller import GenotypeHaplotypeCallerGVCF
from .mutect2 import CallVariantsWithMutect2
from .ref import (CreateCnvBlackListBED, CreateExclusionIntervalListBED,
                  FetchReferenceFASTA, FetchResourceFASTA, FetchResourceFile)
from .samtools import samtools_faidx, samtools_index


class PrepareCanvasGenomeFiles(ShellTask):
    ref_fa_path = luigi.Parameter()
    genomesize_xml_path = luigi.Parameter()
    kmer_fa_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 60

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
        input_genome_fa_path = self.input()[0][0].path
        run_id = Path(input_genome_fa_path).stem
        self.print_log(f'Prepare Canvas genome files:\t{run_id}')
        samtools = self.cf['samtools']
        symlinks = {
            '../{}'.format(Path(p).name): o.path
            for p, o in zip(
                [self.input()[1].path, *[i.path for i in self.input()[2]]],
                self.output()[2:]
            )
        }
        output_genome_fa_path = self.output()[0].path
        output_kmer_fai_path = self.output()[4].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=samtools, cwd=Path(list(symlinks.values())[0]).parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet'], env={'REF_CACHE': '.ref_cache'}
        )
        for s, d in symlinks.items():
            self.run_shell(args=f'ln -s {s} {d}', output_files_or_dirs=d)
        self.run_shell(
            args=(
                f'set -eo pipefail && cut -f 1 {output_kmer_fai_path}'
                + ' | tr "\\n" " "'
                + f' | xargs {samtools} faidx {input_genome_fa_path}'
                + f' > {output_genome_fa_path}'
            ),
            input_files_or_dirs=[
                output_kmer_fai_path, input_genome_fa_path
            ],
            output_files_or_dirs=output_genome_fa_path
        )
        samtools_faidx(
            shelltask=self, samtools=samtools, fa_path=output_genome_fa_path
        )


@requires(CreateCnvBlackListBED, CreateExclusionIntervalListBED)
class CreateCanvasFilterBED(ShellTask):
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath(
                '.'.join([
                    Path(Path(i[0].path).stem).stem for i in self.input()
                ]) + '.bed'
            )
        )

    def run(self):
        output_bed_path = self.output().path
        run_id = Path(output_bed_path).stem
        self.print_log(f'Create a Canvas filter BED:\t{run_id}')
        bedtools = self.cf['bedtools']
        bgzip = self.cf['bgzip']
        input_bed_paths = [i[0].path for i in self.input()]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=bgzip, cwd=self.cf['ref_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                'set -eo pipefail && cat'
                + ''.join([f' <({bgzip} -dc {p})' for p in input_bed_paths])
                + f' | {bedtools} sort'
                + f' | {bedtools} merge > {output_bed_path}'
            ),
            input_files_or_dirs=input_bed_paths,
            output_files_or_dirs=output_bed_path
        )


class CreateUniqueRegionManifest(ShellTask):
    exome_manifest_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 30

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
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
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


class CreateCanvasBAM(ShellTask):
    input_cram_path = luigi.Parameter()
    cram_fa_path = luigi.Parameter()
    kmer_fai_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        bam = Path(self.cf['somatic_cnv_canvas_dir_path']).joinpath(
            Path(self.input_cram_path).stem + '.canvas.bam'
        )
        return [luigi.LocalTarget(f'{bam}{s}') for s in ['', '.bai']]

    def run(self):
        run_id = Path(self.input_cram_path).stem
        self.print_log(f'Create Canvas BAM:\t{run_id}')
        samtools = self.cf['samtools']
        n_cpu = self.cf['n_cpu_per_worker']
        output_bam_path = self.output()[0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=samtools, cwd=Path(output_bam_path).parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet'], env={'REF_CACHE': '.ref_cache'}
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && cut -f 1 {self.kmer_fai_path}'
                + ' | tr "\\n" " "'
                + f' | xargs {samtools} view -@ {n_cpu} -T {self.cram_fa_path}'
                + f' -bS -o {output_bam_path} {self.input_cram_path}'
            ),
            input_files_or_dirs=[
                self.input_cram_path, self.cram_fa_path,
                f'{self.cram_fa_path}.fai', self.kmer_fai_path
            ],
            output_files_or_dirs=output_bam_path
        )
        samtools_index(
            shelltask=self, samtools=samtools, sam_path=output_bam_path,
            n_cpu=n_cpu
        )


@requires(PrepareCRAMTumor, PrepareCRAMNormal, FetchReferenceFASTA,
          PrepareCanvasGenomeFiles, CreateCanvasFilterBED,
          GenotypeHaplotypeCallerGVCF, CallVariantsWithMutect2)
class CallSomaticCopyNumberVariantsWithCanvas(ShellTask):
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    exome_manifest_path = luigi.Parameter(default='')
    priority = 30

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['somatic_cnv_canvas_dir_path']).joinpath(
                    create_matched_id(*[i[0].path for i in self.input()[0:2]])
                )
            )
        ]

    def run(self):
        targets = yield [
            CreateCanvasBAM(
                input_cram_path=self.input()[0][0].path,
                cram_fa_path=self.input()[2][0].path,
                kmer_fai_path=self.input()[3][4].path, cf=self.cf
            ),
            *(
                [
                    CreateUniqueRegionManifest(
                        exome_manifest_path=self.exome_manifest_path,
                        cf=self.cf
                    )
                ] if self.cf['exome'] else list()
            )
        ]
        output_dir_path = self.output()[0].path
        run_id = Path(output_dir_path).name
        self.print_log(f'Call somatic CNVs with Canvas:\t{run_id}')
        canvas = self.cf['Canvas']
        sample_name = self.sample_names[0]
        kmer_fa_path = self.input()[3][3].path
        canvas_genome_dir_path = str(Path(self.input()[3][0].path).parent)
        filter_bed_path = self.input()[4].path
        germline_b_allele_vcf_path = self.input()[5][0].path
        somatic_vcf_path = self.input()[6][0].path
        tumor_bam_path = targets[0][0].path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=canvas, cwd=self.cf['somatic_cnv_canvas_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        if self.cf['exome']:
            manifest_path = targets[1].path
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
        self.run_shell(
            args=f'rm -f {tumor_bam_path} {tumor_bam_path}.bai',
            input_files_or_dirs=tumor_bam_path
        )


class FlagUniqueKmers(ShellTask):
    input_fa_path = luigi.Parameter()
    output_fa_path = luigi.Parameter()
    flaguniquekmers = luigi.Parameter(default='FlagUniqueKmers')

    def output(self):
        return luigi.LocalTarget(self.output_fa_path)

    def run(self):
        run_id = Path(self.input_fa_path).stem
        self.print_log(f'Flag unique k-mers:\t{run_id}')
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            cwd=self.cf['ref_dir_path'], quiet=False
        )
        self.run_shell(
            args=(
                f'set -e && {self.flaguniquekmers}'
                + f' {self.input_fa_path} {self.output_fa_path}'
            ),
            input_files_or_dirs=self.input_fa_path,
            output_files_or_dirs=self.output_fa_path
        )


if __name__ == '__main__':
    luigi.run()
