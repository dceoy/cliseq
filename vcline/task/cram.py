#!/usr/bin/env python

import sys
from pathlib import Path

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
        if self.cram_list:
            input_sam = Path(self.cram_list[0])
            return SamtoolsView(
                input_sam_path=str(input_sam),
                output_sam_path=str(
                    input_sam.parent.joinpath(f'{input_sam.stem}.cram')
                ),
                fa_path=self.ref_fa_path, samtools=self.cf['samtools'],
                n_cpu=self.n_cpu, remove_input=False,
                index_sam=True, sh_config=self.sh_config
            )
        else:
            return PrepareAnalysisReadyCram(
                fq_paths=self.fq_list[0], read_group=self.read_groups[0],
                sample_name=self.sample_names[0], ref_fa_path=self.ref_fa_path,
                known_sites_vcf_paths=[
                    self.dbsnp_vcf_path, self.mills_indel_vcf_path,
                    self.known_indel_vcf_path,
                ],
                cf={
                    'metrics_collectors': list(),
                    **{k: v for k, v in self.cf if k != 'metrics_collectors'}
                },
                n_cpu=self.n_cpu, memory_mb=self.memory_mb
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
        if self.cram_list:
            input_sam = Path(self.cram_list[1])
            return SamtoolsView(
                input_sam_path=str(input_sam),
                output_sam_path=str(
                    input_sam.parent.joinpath(f'{input_sam.stem}.cram')
                ),
                fa_path=self.ref_fa_path, samtools=self.cf['samtools'],
                n_cpu=self.n_cpu, remove_input=False,
                index_sam=True, sh_config=self.sh_config
            )
        else:
            return PrepareAnalysisReadyCram(
                fq_paths=self.fq_list[1], read_group=self.read_groups[1],
                sample_name=self.sample_names[1], ref_fa_path=self.ref_fa_path,
                known_sites_vcf_paths=[
                    self.dbsnp_vcf_path, self.mills_indel_vcf_path,
                    self.known_indel_vcf_path,
                ],
                cf={
                    'metrics_collectors': list(),
                    **{k: v for k, v in self.cf if k != 'metrics_collectors'}
                },
                n_cpu=self.n_cpu, memory_mb=self.memory_mb
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
