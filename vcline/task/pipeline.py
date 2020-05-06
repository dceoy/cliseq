#!/usr/bin/env python

import logging
import os
import sys
from itertools import chain
from pathlib import Path

import luigi
from luigi.tools import deps_tree

from ..cli.util import create_matched_id
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .base import ShellTask
from .callcopyratiosegments import CallCopyRatioSegmentsMatched
from .canvas import CallSomaticCopyNumberVariantsWithCanvas
from .delly import CallStructualVariantsWithDelly
from .funcotator import FuncotateSegments, FuncotateVariants
from .haplotypecaller import FilterVariantTranches
from .manta import CallStructualVariantsWithManta
from .msisensor import ScoreMSIWithMSIsensor
from .mutect2 import FilterMutectCalls
from .snpeff import AnnotateVariantsWithSnpEff
from .strelka import (CallGermlineVariantsWithStrelka,
                      CallSomaticVariantsWithStrelka)


class RunVariantCaller(luigi.Task):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    hapmap_vcf_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    cnv_black_list_path = luigi.Parameter()
    genomesize_xml_path = luigi.Parameter()
    kmer_fa_path = luigi.Parameter()
    exome_manifest_path = luigi.Parameter(default='')
    funcotator_somatic_tar_path = luigi.Parameter(default='')
    funcotator_germline_tar_path = luigi.Parameter(default='')
    snpeff_config_path = luigi.Parameter(default='')
    cf = luigi.DictParameter()
    caller = luigi.Parameter()
    annotators = luigi.ListParameter(default=list())
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = luigi.IntParameter(default=1000)

    def requires(self):
        if 'germline_short_variant.gatk' == self.caller:
            return FilterVariantTranches(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                hapmap_vcf_path=self.hapmap_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_short_variant.gatk' == self.caller:
            return FilterMutectCalls(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                gnomad_vcf_path=self.gnomad_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_structual_variant.manta' == self.caller:
            return CallStructualVariantsWithManta(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_short_variant.strelka' == self.caller:
            return CallSomaticVariantsWithStrelka(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'germline_short_variant.strelka' == self.caller:
            return CallGermlineVariantsWithStrelka(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_structual_variant.delly' == self.caller:
            return CallStructualVariantsWithDelly(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_copy_number_variation.gatk' == self.caller:
            return CallCopyRatioSegmentsMatched(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cnv_black_list_path=self.cnv_black_list_path, cf=self.cf
            )
        elif 'somatic_copy_number_variation.canvas' == self.caller:
            return CallSomaticCopyNumberVariantsWithCanvas(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                gnomad_vcf_path=self.gnomad_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cnv_black_list_path=self.cnv_black_list_path,
                genomesize_xml_path=self.genomesize_xml_path,
                kmer_fa_path=self.kmer_fa_path,
                exome_manifest_path=self.exome_manifest_path, cf=self.cf
            )
        elif 'somatic_msi.msisensor' == self.caller:
            return ScoreMSIWithMSIsensor(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        else:
            raise ValueError(f'invalid caller: {self.caller}')

    def output(self):
        output_pos = list(
            chain.from_iterable([
                [
                    Path(self.cf[f'postproc_{k}_dir_path']).joinpath(
                        (Path(p).name + f'.{k}.tsv') if p.endswith('.seg')
                        else (Path(Path(p).stem).stem + f'.norm.{k}.vcf.gz')
                    ) for p in v
                ] for k, v in self._find_annotation_targets().items()
            ])
        )
        return (
            [luigi.LocalTarget(o) for o in output_pos]
            if output_pos else self.input()
        )

    def _find_annotation_targets(self):
        input_paths = [i.path for i in self.input()]
        if 'somatic_structual_variant.delly' == self.caller:
            suffix_dict = {'funcotator': None, 'snpeff': '.vcf.gz'}
        elif 'somatic_structual_variant.manta' == self.caller:
            suffix_dict = {
                'funcotator': '.manta.somaticSV.vcf.gz', 'snpeff': '.vcf.gz'
            }
        else:
            suffix_dict = {
                'funcotator': ('.vcf.gz', '.called.seg'), 'snpeff': '.vcf.gz'
            }
        return {
            k: (
                [p for p in input_paths if v and p.endswith(v)]
                if k in self.annotators else list()
            ) for k, v in suffix_dict.items()
        }

    def run(self):
        logger = logging.getLogger(__name__)
        logger.debug('Task tree:' + os.linesep + deps_tree.print_tree(self))
        funcotator_common_kwargs = {
            'data_src_tar_path': (
                self.funcotator_germline_tar_path
                if self.caller.startswith('germline_')
                else self.funcotator_somatic_tar_path
            ),
            'ref_fa_path': self.ref_fa_path, 'cf': self.cf
        }
        for k, v in self._find_annotation_targets().items():
            if k == 'funcotator':
                yield [
                    (
                        FuncotateSegments(
                            input_seg_path=p, **funcotator_common_kwargs
                        ) if p.endswith('.seg') else FuncotateVariants(
                            input_vcf_path=p, normalize_vcf=self.normalize_vcf,
                            **funcotator_common_kwargs
                        )
                    ) for p in v
                ]
            elif k == 'snpeff':
                yield [
                    AnnotateVariantsWithSnpEff(
                        input_vcf_path=p,
                        snpeff_config_path=self.snpeff_config_path,
                        ref_fa_path=self.ref_fa_path, cf=self.cf,
                        normalize_vcf=self.normalize_vcf
                    ) for p in v
                ]


class PrintEnvVersions(ShellTask):
    log_dir_path = luigi.Parameter()
    command_paths = luigi.ListParameter(default=list())
    run_id = luigi.Parameter(default='env')
    priority = luigi.IntParameter(default=sys.maxsize)
    __is_completed = False

    def complete(self):
        return self.__is_completed

    def run(self):
        python = sys.executable
        self.print_log(f'Print environment versions: {python}')
        self.setup_shell(
            run_id=self.run_id, log_dir_path=self.log_dir_path,
            commands=[python, *self.command_paths], quiet=False
        )
        self.run_shell(
            args=[
                f'{python} -m pip --version',
                f'{python} -m pip freeze --no-cache-dir'
            ]
        )
        self.__is_completed = True


class RemoveBAMs(ShellTask):
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    ref_fa_path = luigi.Parameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    hapmap_vcf_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    cnv_black_list_path = luigi.Parameter()
    genomesize_xml_path = luigi.Parameter()
    kmer_fa_path = luigi.Parameter()
    exome_manifest_path = luigi.Parameter(default='')
    funcotator_somatic_tar_path = luigi.Parameter()
    funcotator_germline_tar_path = luigi.Parameter()
    snpeff_config_path = luigi.Parameter()
    cf = luigi.DictParameter()
    callers = luigi.ListParameter()
    annotators = luigi.ListParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = luigi.IntParameter(default=0)
    __is_completed = False

    def requires(self):
        kwargs = {
            k: getattr(self, k) for k in [
                'fq_list', 'read_groups', 'cram_list', 'sample_names',
                'ref_fa_path', 'dbsnp_vcf_path', 'mills_indel_vcf_path',
                'known_indel_vcf_path', 'cf'
            ]
        }
        return [PrepareCRAMTumor(**kwargs), PrepareCRAMNormal(**kwargs)]

    def complete(self):
        return self.__is_completed

    def run(self):
        align_dir = Path(self.cf['align_dir_path'])
        bams = [
            align_dir.joinpath(Path(i[0].path).stem + '.bam')
            for i in self.input()
        ]
        if any([b.is_file() for b in bams]):
            run_id = create_matched_id(*[b.name for b in bams])
            self.print_log(f'Remove BAM files:\t{run_id}')
            self.setup_shell(
                run_id=run_id, log_dir_path=self.cf['log_dir_path'],
                cwd=self.cf['align_dir_path'],
                remove_if_failed=self.cf['remove_if_failed']
            )
            self.run_shell(
                args=('rm -f' + ''.join([f' {b} {b}.bai' for b in bams]))
            )
        if not list(align_dir.iterdir()):
            align_dir.rmdir()
        self.__is_completed = True


if __name__ == '__main__':
    luigi.run()
