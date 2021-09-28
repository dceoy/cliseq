#!/usr/bin/env python

import logging
import os
from pathlib import Path
from socket import gethostname

import luigi
from luigi.tools import deps_tree
from vanqc.task.bcftools import CollectVcfStats, NormalizeVcf
from vanqc.task.gatk import (AnnotateSegWithFuncotateSegments,
                             AnnotateVariantsWithFuncotator)
from vanqc.task.picard import CollectVariantCallingMetrics
from vanqc.task.snpeff import AnnotateVariantsWithSnpeff
from vanqc.task.vep import AnnotateVariantsWithEnsemblVep

from .callcopyratiosegments import CallCopyRatioSegmentsMatched
from .core import VclineTask
from .delly import CallSomaticStructualVariantsWithDelly
from .haplotypecaller import FilterVariantTranches
from .manta import CallSomaticStructualVariantsWithManta
from .msisensor import ScoreMsiWithMsisensor
from .mutect2 import FilterMutectCalls
from .strelka import (CallGermlineVariantsWithStrelka,
                      CallSomaticVariantsWithStrelka)


class PrintEnvVersions(VclineTask):
    command_paths = luigi.ListParameter(default=list())
    run_id = luigi.Parameter(default=gethostname())
    sh_config = luigi.DictParameter(default=dict())
    __is_completed = False

    def complete(self):
        return self.__is_completed

    def run(self):
        self.print_log(f'Print environment versions:\t{self.run_id}')
        self.setup_shell(
            run_id=self.run_id, commands=self.command_paths, **self.sh_config
        )
        self.print_env_versions()
        self.__is_completed = True


class RunVariantCaller(luigi.Task):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    cf = luigi.DictParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    hapmap_vcf_path = luigi.Parameter(default='')
    kg_snps_vcf_path = luigi.Parameter(default='')
    gnomad_vcf_path = luigi.Parameter(default='')
    cnv_blacklist_path = luigi.Parameter(default='')
    funcotator_somatic_data_dir_path = luigi.Parameter(default='')
    funcotator_germline_data_dir_path = luigi.Parameter(default='')
    snpeff_db_data_dir_path = luigi.Parameter(default='')
    vep_cache_data_dir_path = luigi.Parameter(default='')
    caller = luigi.Parameter(default='somatic_snv_indel.gatk')
    metrics_collectors = luigi.ListParameter(default=list())
    annotators = luigi.ListParameter(default=list())
    normalize_vcf = luigi.BoolParameter(default=True)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = luigi.IntParameter(default=1000)

    def requires(self):
        if 'somatic_snv_indel.gatk' == self.caller:
            assert bool(self.gnomad_vcf_path)
            return FilterMutectCalls(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                gnomad_vcf_path=self.gnomad_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )
        elif 'germline_snv_indel.gatk' == self.caller:
            assert bool(self.hapmap_vcf_path and self.kg_snps_vcf_path)
            return FilterVariantTranches(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                hapmap_vcf_path=self.hapmap_vcf_path,
                kg_snps_vcf_path=self.kg_snps_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )
        elif 'somatic_sv.manta' == self.caller:
            return CallSomaticStructualVariantsWithManta(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
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
                cf=self.cf, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
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
                cf=self.cf, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )
        elif 'somatic_sv.delly' == self.caller:
            return CallSomaticStructualVariantsWithDelly(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )
        elif 'somatic_cnv.gatk' == self.caller:
            assert bool(self.cnv_blacklist_path)
            return CallCopyRatioSegmentsMatched(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                hapmap_vcf_path=self.hapmap_vcf_path,
                kg_snps_vcf_path=self.kg_snps_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cnv_blacklist_path=self.cnv_blacklist_path, cf=self.cf,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )
        elif 'somatic_msi.msisensor' == self.caller:
            return ScoreMsiWithMsisensor(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )
        else:
            raise ValueError(f'invalid caller: {self.caller}')

    def output(self):
        postproc_dir = Path(self.cf['postproc_dir_path'])
        norm_tag = ('.norm' if self.normalize_vcf else '')
        return (
            [
                luigi.LocalTarget(
                    postproc_dir.joinpath(a).joinpath(
                        (
                            Path(Path(p).stem).stem + {
                                'norm': f'{norm_tag}.vcf.gz',
                                'funcotator': f'{norm_tag}.funcotator.vcf.gz',
                                'snpeff': f'{norm_tag}.snpeff.vcf.gz',
                                'vep': f'{norm_tag}.vep.txt.gz',
                                'bcftools': '.vcf.stats.txt',
                                'picard': '.CollectVariantCallingMetrics.txt'
                            }[a]
                        ) if p.endswith('.vcf.gz')
                        else (Path(p).stem + f'.{a}.seg.tsv')
                    )
                ) for p, a in self._generate_postproc_targets()
            ] or self.input()
        )

    def _generate_postproc_targets(self):
        for i in self.input():
            p = i.path
            if p.endswith('.vcf.gz'):
                for c in self.metrics_collectors:
                    yield p, c
                if self.normalize_vcf and not self.annotators:
                    yield p, 'norm'
                else:
                    for a in self.annotators:
                        if (a == 'funcotator'
                                and self.caller.startswith('somatic_sv.')):
                            pass
                        else:
                            yield p, a
            elif p.endswith('.called.seg') and 'funcotator' in self.annotators:
                yield p, 'funcotator'

    def run(self):
        postproc_dir = Path(self.cf['postproc_dir_path'])
        norm_dir = postproc_dir.joinpath('norm')
        for p, a in self._generate_postproc_targets():
            if a == 'bcftools':
                yield CollectVcfStats(
                    input_vcf_path=p, fa_path=self.ref_fa_path,
                    dest_dir_path=str(postproc_dir.joinpath(a)),
                    bcftools=self.cf['bcftools'],
                    plot_vcfstats=self.cf['plot_vcfstats'], n_cpu=self.n_cpu,
                    sh_config=self.sh_config
                )
            elif a == 'picard':
                yield CollectVariantCallingMetrics(
                    input_vcf_path=p, fa_path=self.ref_fa_path,
                    dbsnp_vcf_path=self.dbsnp_vcf_path,
                    dest_dir_path=str(postproc_dir.joinpath(a)),
                    picard=self.cf['gatk'], n_cpu=self.n_cpu,
                    memory_mb=self.memory_mb, sh_config=self.sh_config
                )
            elif a == 'funcotator':
                data_src_dir_path = (
                    self.funcotator_germline_data_dir_path
                    if self.caller.startswith('germline_')
                    else self.funcotator_somatic_data_dir_path
                )
                assert bool(data_src_dir_path)
                if p.endswith('.seg'):
                    yield AnnotateSegWithFuncotateSegments(
                        input_seg_path=p, fa_path=self.ref_fa_path,
                        data_src_dir_path=data_src_dir_path,
                        ref_version=self.cf['ucsc_hg_version'],
                        dest_dir_path=str(postproc_dir.joinpath(a)),
                        gatk=self.cf['gatk'], n_cpu=self.n_cpu,
                        memory_mb=self.memory_mb, sh_config=self.sh_config
                    )
                else:
                    yield AnnotateVariantsWithFuncotator(
                        input_vcf_path=p, fa_path=self.ref_fa_path,
                        data_src_dir_path=data_src_dir_path,
                        ref_version=self.cf['ucsc_hg_version'],
                        dest_dir_path=str(postproc_dir.joinpath(a)),
                        normalize_vcf=self.normalize_vcf,
                        norm_dir_path=str(norm_dir),
                        bcftools=self.cf['bcftools'], gatk=self.cf['gatk'],
                        n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                        sh_config=self.sh_config
                    )
            elif a == 'snpeff':
                assert bool(self.snpeff_db_data_dir_path)
                yield AnnotateVariantsWithSnpeff(
                    input_vcf_path=p, fa_path=self.ref_fa_path,
                    db_data_dir_path=self.snpeff_db_data_dir_path,
                    dest_dir_path=str(postproc_dir.joinpath(a)),
                    normalize_vcf=self.normalize_vcf,
                    norm_dir_path=str(norm_dir), bcftools=self.cf['bcftools'],
                    snpeff=self.cf['snpEff'], bgzip=self.cf['bgzip'],
                    tabix=self.cf['tabix'], n_cpu=self.n_cpu,
                    memory_mb=self.memory_mb, sh_config=self.sh_config
                )
            elif a == 'vep':
                assert bool(self.vep_cache_data_dir_path)
                yield AnnotateVariantsWithEnsemblVep(
                    input_vcf_path=p, fa_path=self.ref_fa_path,
                    cache_data_dir_path=self.vep_cache_data_dir_path,
                    dest_dir_path=str(postproc_dir.joinpath(a)),
                    normalize_vcf=self.normalize_vcf,
                    norm_dir_path=str(norm_dir), bcftools=self.cf['bcftools'],
                    vep=self.cf['vep'], pigz=self.cf['pigz'],
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                    sh_config=self.sh_config
                )
            elif a == 'norm' and not self.annotators:
                yield NormalizeVcf(
                    input_vcf_path=p, fa_path=self.ref_fa_path,
                    dest_dir_path=str(norm_dir), bcftools=self.cf['bcftools'],
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                    sh_config=self.sh_config
                )
        logger = logging.getLogger(__name__)
        logger.debug('Task tree:' + os.linesep + deps_tree.print_tree(self))


if __name__ == '__main__':
    luigi.run()
