#!/usr/bin/env python

import sys

import luigi
from ftarc.task.controller import PrepareAnalysisReadyCram
from ftarc.task.samtools import SamtoolsView
from luigi.util import requires


class PrepareCramTumor(luigi.WrapperTask):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = luigi.IntParameter(default=100)

    def requires(self):
        return (
            SamtoolsView(
                input_sam_path=self.cram_list[0],
                output_sam_path=self.output()[0].path,
                fa_path=self.ref_fa_path, samtools=self.cf['samtools'],
                n_cpu=self.n_cpu, remove_input=False,
                index_sam=True, sh_config=self.sh_config
            ) if self.cram_list
            else PrepareAnalysisReadyCram(
                fq_paths=self.fq_list[0], read_group=self.read_groups[0],
                sample_name=self.sample_names[0], ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path, cf=self.cf,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb
            )
        )

    def output(self):
        return self.input()


class PrepareCramNormal(luigi.WrapperTask):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = luigi.IntParameter(default=100)

    def requires(self):
        return (
            SamtoolsView(
                input_sam_path=self.cram_list[1],
                output_sam_path=self.output()[0].path,
                fa_path=self.ref_fa_path, samtools=self.cf['samtools'],
                n_cpu=self.n_cpu, remove_input=False,
                index_sam=True, sh_config=self.sh_config
            ) if self.cram_list
            else PrepareAnalysisReadyCram(
                fq_paths=self.fq_list[1], read_group=self.read_groups[1],
                sample_name=self.sample_names[1], ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path, cf=self.cf,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb
            )
        )

    def output(self):
        return self.input()


@requires(PrepareCramTumor, PrepareCramNormal)
class PrepareCramsMatched(luigi.WrapperTask):
    priority = luigi.IntParameter(default=sys.maxsize)

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
