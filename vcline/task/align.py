#!/usr/bin/env python

from pathlib import Path

import luigi
from ftarc.task.controller import PrepareAnalysisReadyCRAM
from ftarc.task.resource import FetchReferenceFASTA
from ftarc.task.samtools import SamtoolsIndex, SamtoolsView


class PrepareCRAMTumor(luigi.Task):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = luigi.IntParameter(default=100)

    def requires(self):
        if not self.cram_list:
            return PrepareAnalysisReadyCRAM(
                fq_paths=self.fq_list[0], read_group=self.read_groups[0],
                sample_name=self.sample_names[0], ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path, cf=self.cf
            )
        else:
            return super().requires()

    def output(self):
        if not self.cram_list:
            return self.input()
        else:
            cram_path = (
                self.cram_list[0] if self.cram_list[0].endswith('.cram')
                else str(
                    Path(self.cf['align_dir_path']).joinpath(
                        self.sample_names[0]
                    ).joinpath(Path(self.cram_list[0]).stem + '.cram')
                )
            )
            return [luigi.LocalTarget(cram_path + s) for s in ['', '.crai']]

    def run(self):
        if not self.cram_list:
            pass
        elif self.cram_list[0].endswith('.cram'):
            yield SamtoolsIndex(
                sam_path=self.cram_list[0], samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'],
                log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed'],
                quiet=self.cf['quiet']
            )
        else:
            ref_target = yield FetchReferenceFASTA(
                ref_fa_path=self.ref_fa_path, cf=self.cf
            )
            yield SamtoolsView(
                input_sam_path=self.cram_list[0],
                output_sam_path=self.output()[0].path,
                fa_path=ref_target[0].path, samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'], remove_input=False,
                index_sam=True, log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed'],
                quiet=self.cf['quiet']
            )


class PrepareCRAMNormal(luigi.Task):
    ref_fa_path = luigi.Parameter()
    fq_list = luigi.ListParameter()
    cram_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    cf = luigi.DictParameter()
    priority = luigi.IntParameter(default=100)

    def requires(self):
        if not self.cram_list:
            return PrepareAnalysisReadyCRAM(
                fq_paths=self.fq_list[1], read_group=self.read_groups[1],
                sample_name=self.sample_names[1], ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path, cf=self.cf
            )
        else:
            return super().requires()

    def output(self):
        if not self.cram_list:
            return self.input()
        else:
            cram_path = (
                self.cram_list[1] if self.cram_list[1].endswith('.cram')
                else str(
                    Path(self.cf['align_dir_path']).joinpath(
                        self.sample_names[1]
                    ).joinpath(Path(self.cram_list[1]).stem + '.cram')
                )
            )
            return [luigi.LocalTarget(cram_path + s) for s in ['', '.crai']]

    def run(self):
        if not self.cram_list:
            pass
        elif self.cram_list[1].endswith('.cram'):
            yield SamtoolsIndex(
                sam_path=self.cram_list[1], samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'],
                log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed'],
                quiet=self.cf['quiet']
            )
        else:
            ref_target = yield FetchReferenceFASTA(
                ref_fa_path=self.ref_fa_path, cf=self.cf
            )
            yield SamtoolsView(
                input_sam_path=self.cram_list[1],
                output_sam_path=self.output()[0].path,
                fa_path=ref_target[0].path, samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'], remove_input=False,
                index_sam=True, log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed'],
                quiet=self.cf['quiet']
            )


if __name__ == '__main__':
    luigi.run()
