#!/usr/bin/env python

import luigi

from .align import ApplyBQSR
from .ref import CreateFASTAIndex, FetchGenomeFASTA


class CallVariants(luigi.WrapperTask):
    ref_fa_list = luigi.ListParameter()
    fq_dict = luigi.DictParameter()
    cf = luigi.DictParameter()
    priority = 10

    def requires(self):
        return [
            FetchGenomeFASTA(ref_fa_list=self.ref_fa_list, cf=self.cf),
            CreateFASTAIndex(ref_fa_list=self.ref_fa_list, cf=self.cf),
            *[
                ApplyBQSR(fq_paths=v, ref_fa_list=self.ref_fa_list, cf=self.cf)
                for v in self.fq_dict.values()
            ]
        ]
