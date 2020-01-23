#!/usr/bin/env python

import luigi

from .align import ApplyBQSR


class PrepareCRAMs(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_dict = luigi.DictParameter()
    known_site_vcf_paths = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 10

    def requires(self):
        return [
            ApplyBQSR(
                fq_paths=v, ref_fa_paths=self.ref_fa_paths,
                known_site_vcf_paths=self.known_site_vcf_paths, cf=self.cf
            ) for v in self.fq_dict.values()
        ]
