#!/usr/bin/env python

import logging
import os
import sys
from itertools import chain
from pathlib import Path

import luigi
from luigi.tools import deps_tree
from luigi.util import requires

from .base import ShellTask
from .bcftools import NormalizeVCF
from .callcopyratiosegments import CallCopyRatioSegmentsMatched
from .delly import CallStructualVariantsWithDelly
from .funcotator import FuncotateSegments, FuncotateVariants
from .haplotypecaller import FilterVariantTranches
from .manta import CallStructualVariantsWithManta
from .msisensor import ScoreMSIWithMSIsensor
from .mutect2 import FilterMutectCalls
from .ref import (CreateBWAIndices, CreateCnvBlackListBED,
                  CreateEvaluationIntervalListBED,
                  CreateExclusionIntervalListBED, CreateGnomadBiallelicSnpVCF,
                  CreateSequenceDictionary, FetchCnvBlackList, FetchDbsnpVCF,
                  FetchEvaluationIntervalList, FetchGnomadVCF, FetchHapmapVCF,
                  FetchKnownIndelVCF, FetchMillsIndelVCF, FetchReferenceFASTA,
                  PreprocessIntervals, ScanMicrosatellites,
                  UncompressEvaluationIntervalListBED)
from .snpeff import AnnotateVariantsWithSnpEff
from .strelka import (CallGermlineVariantsWithStrelka,
                      CallSomaticVariantsWithStrelka)


class RunVariantCaller(luigi.Task):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter(default='')
    mills_indel_vcf_path = luigi.Parameter(default='')
    known_indel_vcf_path = luigi.Parameter(default='')
    hapmap_vcf_path = luigi.Parameter(default='')
    gnomad_vcf_path = luigi.Parameter(default='')
    evaluation_interval_path = luigi.Parameter(default='')
    cnv_blacklist_path = luigi.Parameter(default='')
    funcotator_somatic_tar_path = luigi.Parameter(default='')
    funcotator_germline_tar_path = luigi.Parameter(default='')
    snpeff_config_path = luigi.Parameter(default='')
    cf = luigi.DictParameter()
    caller = luigi.Parameter()
    annotators = luigi.ListParameter(default=list())
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = luigi.IntParameter(default=1000)

    def requires(self):
        if 'germline_snv_indel.gatk' == self.caller:
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
        elif 'somatic_snv_indel.gatk' == self.caller:
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
        elif 'somatic_sv.manta' == self.caller:
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
        elif 'somatic_snv_indel.strelka' == self.caller:
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
        elif 'germline_snv_indel.strelka' == self.caller:
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
        elif 'somatic_sv.delly' == self.caller:
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
        elif 'somatic_cnv.gatk' == self.caller:
            return CallCopyRatioSegmentsMatched(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cnv_blacklist_path=self.cnv_blacklist_path, cf=self.cf
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
        output_files = list(
            chain.from_iterable([
                [
                    Path(self.cf['postproc_dir_path']).joinpath(
                        'normalization' if a == 'normalize' else 'annotation'
                    ).joinpath(
                        (Path(p).name + f'.{a}.tsv') if p.endswith('.seg')
                        else (
                            Path(Path(p).stem).stem
                            + ('.norm' if self.normalize_vcf else '')
                            + ('' if a == 'normalize' else f'.{a}')
                            + '.vcf.gz'
                        )
                    )
                ] for p, a in self._generate_annotation_targets()
            ])
        )
        return (
            [luigi.LocalTarget(o) for o in output_files]
            if output_files else self.input()
        )

    def _generate_annotation_targets(self):
        for i in self.input():
            p = i.path
            if p.endswith('.called.seg') and 'funcotator' in self.annotators:
                yield p, 'funcotator'
            elif p.endswith('.vcf.gz'):
                if self.normalize_vcf:
                    yield p, 'normalize'
                if ('funcotator' in self.annotators
                        and not self.caller.startswith('somatic_sv.')):
                    yield p, 'funcotator'
                if 'snpeff' in self.annotators:
                    yield p, 'snpeff'

    def run(self):
        for p, a in self._generate_annotation_targets():
            if a == 'funcotator':
                data_src_tar_path = (
                    self.funcotator_germline_tar_path
                    if self.caller.startswith('germline_')
                    else self.funcotator_somatic_tar_path
                )
                if p.endswith('.seg'):
                    yield FuncotateSegments(
                        input_seg_path=p, ref_fa_path=self.ref_fa_path,
                        data_src_tar_path=data_src_tar_path, cf=self.cf
                    )
                else:
                    yield FuncotateVariants(
                        input_vcf_path=p, ref_fa_path=self.ref_fa_path,
                        data_src_tar_path=data_src_tar_path, cf=self.cf,
                        normalize_vcf=self.normalize_vcf
                    )
            elif a == 'snpeff':
                yield AnnotateVariantsWithSnpEff(
                    input_vcf_path=p,
                    snpeff_config_path=self.snpeff_config_path,
                    ref_fa_path=self.ref_fa_path, cf=self.cf,
                    normalize_vcf=self.normalize_vcf
                )
            elif a == 'normalize' and not self.annotators:
                ref_target = yield FetchReferenceFASTA(
                    ref_fa_path=self.ref_fa_path, cf=self.cf
                )
                yield NormalizeVCF(
                    input_vcf_path=p, fa_path=ref_target[0].path,
                    dest_dir_path=str(Path(self.output()[0].path).parent),
                    n_cpu=self.cf['n_cpu_per_worker'],
                    memory_mb=self.cf['memory_mb_per_worker'],
                    bcftools=self.cf['bcftools'],
                    log_dir_path=self.cf['log_dir_path'],
                    remove_if_failed=self.cf['remove_if_failed'],
                    quiet=self.cf['quiet']
                )
        logger = logging.getLogger(__name__)
        logger.debug('Task tree:' + os.linesep + deps_tree.print_tree(self))


class PrintEnvVersions(ShellTask):
    log_dir_path = luigi.Parameter()
    command_paths = luigi.ListParameter(default=list())
    run_id = luigi.Parameter(default='env')
    quiet = luigi.BoolParameter(default=False)
    priority = luigi.IntParameter(default=sys.maxsize)
    __is_completed = False

    def complete(self):
        return self.__is_completed

    def run(self):
        python = sys.executable
        self.print_log(f'Print environment versions: {python}')
        version_files = [
            Path('/proc/version'),
            *[
                o for o in Path('/etc').iterdir()
                if o.name.endswith(('-release', '_version'))
            ]
        ]
        self.setup_shell(
            run_id=self.run_id, log_dir_path=self.log_dir_path,
            commands=[python, *self.command_paths], quiet=self.quiet
        )
        self.run_shell(
            args=[
                f'{python} -m pip --version',
                f'{python} -m pip freeze --no-cache-dir'
            ]
        )
        self.run_shell(
            args=[
                'uname -a',
                *[f'cat {o}' for o in version_files if o.is_file()]
            ]
        )
        self.__is_completed = True


@requires(FetchReferenceFASTA, CreateBWAIndices, CreateSequenceDictionary,
          FetchDbsnpVCF, FetchMillsIndelVCF, FetchKnownIndelVCF,
          FetchEvaluationIntervalList, CreateEvaluationIntervalListBED,
          CreateExclusionIntervalListBED, FetchHapmapVCF, FetchGnomadVCF,
          CreateGnomadBiallelicSnpVCF, FetchCnvBlackList,
          CreateCnvBlackListBED, ScanMicrosatellites,
          UncompressEvaluationIntervalListBED)
class PreprocessResources(luigi.Task):
    ref_fa_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    cnv_blacklist_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = luigi.IntParameter(default=sys.maxsize)
    __is_completed = False

    def complete(self):
        return self.__is_completed

    def run(self):
        yield [
            PreprocessIntervals(
                ref_fa_path=self.ref_fa_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cnv_blacklist_path=self.cnv_blacklist_path,
                cf={
                    k: (bool(i) if k == 'exome' else v)
                    for k, v in self.cf.items()
                }
            ) for i in range(2)
        ]
        self.__is_completed = True


if __name__ == '__main__':
    luigi.run()
