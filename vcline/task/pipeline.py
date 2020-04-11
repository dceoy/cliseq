#!/usr/bin/env python

import logging
import os
from pathlib import Path
from pprint import pformat

import luigi
from luigi.tools import deps_tree

from ..cli.util import fetch_executable, parse_fq_id, read_config_yml
from .align import PrepareCRAMs
from .base import BaseTask
from .callcopyratiosegments import CallCopyRatioSegmentsTumor
from .delly import CallStructualVariantsWithDelly
from .funcotator import AnnotateVariantsWithFuncotator
from .haplotypecaller import FilterVariantTranches
from .lumpy import CallStructualVariantsWithLumpy
from .manta import CallStructualVariantsWithManta
from .mutect2 import FilterMutectCalls
from .strelka import (CallGermlineVariantsWithStrelka,
                      CallSomaticVariantsWithStrelka)


class RunVariantCaller(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    hapmap_vcf_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    cnv_black_list_path = luigi.Parameter()
    funcotator_data_src_tar_path = luigi.Parameter()
    cf = luigi.DictParameter()
    caller_mode = luigi.ListParameter()
    annotators = luigi.ListParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = 10

    def requires(self):
        if 'germline_short_variant.gatk' == self.caller_mode:
            return FilterVariantTranches(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                hapmap_vcf_path=self.hapmap_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_short_variant.gatk' == self.caller_mode:
            return FilterMutectCalls(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                gnomad_vcf_path=self.gnomad_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_structual_variant.manta' == self.caller_mode:
            return CallStructualVariantsWithManta(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_short_variant.strelka' == self.caller_mode:
            return CallSomaticVariantsWithStrelka(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'germline_short_variant.strelka' == self.caller_mode:
            return CallGermlineVariantsWithStrelka(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_structual_variant.delly' == self.caller_mode:
            return CallStructualVariantsWithDelly(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_structual_variant.lumpy' == self.caller_mode:
            return CallStructualVariantsWithLumpy(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_copy_number_variation.gatk' == self.caller_mode:
            return CallCopyRatioSegmentsTumor(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cnv_black_list_path=self.cnv_black_list_path, cf=self.cf
            )
        else:
            raise ValueError(f'invalid caller mode: {self.caller_mode}')

    def output(self):
        vcf_paths = [
            i.path for i in self.input() if (
                i.path.endswith(('.vcf.gz', '.bcf'))
                and self.caller_mode in Path(i.path).name
            )
        ]
        if vcf_paths and 'funcotator' in self.annotators:
            return [
                AnnotateVariantsWithFuncotator(
                    input_vcf_path=p,
                    data_src_tar_path=self.funcotator_data_src_tar_path,
                    ref_fa_paths=self.ref_fa_paths,
                    cf=self.cf, normalize_vcf=self.normalize_vcf
                ).output() for p in vcf_paths
            ]
        else:
            return self.input()


class CallVariants(luigi.WrapperTask):
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    ref_fa_paths = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    mills_indel_vcf_path = luigi.Parameter()
    known_indel_vcf_path = luigi.Parameter()
    hapmap_vcf_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    cnv_black_list_path = luigi.Parameter()
    funcotator_somatic_tar_path = luigi.Parameter()
    funcotator_germline_tar_path = luigi.Parameter()
    cf = luigi.DictParameter()
    caller_modes = luigi.ListParameter()
    annotators = luigi.ListParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = luigi.IntParameter(default=100)

    def requires(self):
        if self.caller_modes:
            return [
                RunVariantCaller(
                    fq_list=self.fq_list, read_groups=self.read_groups,
                    sample_names=self.sample_names,
                    ref_fa_paths=self.ref_fa_paths,
                    dbsnp_vcf_path=self.dbsnp_vcf_path,
                    mills_indel_vcf_path=self.mills_indel_vcf_path,
                    known_indel_vcf_path=self.known_indel_vcf_path,
                    hapmap_vcf_path=self.hapmap_vcf_path,
                    gnomad_vcf_path=self.gnomad_vcf_path,
                    evaluation_interval_path=self.evaluation_interval_path,
                    cnv_black_list_path=self.cnv_black_list_path,
                    funcotator_data_src_tar_path=(
                        self.funcotator_germline_tar_path
                        if m.endswith('.germline') else
                        self.funcotator_somatic_tar_path
                    ),
                    cf=self.cf, caller_mode=m,
                    annotators=self.annotators,
                    normalize_vcf=self.normalize_vcf
                ) for m in self.caller_modes
            ]
        else:
            return PrepareCRAMs(
                fq_list=self.fq_list, read_groups=self.read_groups,
                sample_names=self.sample_names, ref_fa_paths=self.ref_fa_paths,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path, cf=self.cf
            )

    def output(self):
        return self.input()


class RunAnalyticalPipeline(BaseTask):
    config_yml_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    log_dir_path = luigi.Parameter(default='log')
    ref_dir_path = luigi.Parameter(default='')
    n_cpu_per_worker = luigi.IntParameter(default=1)
    memory_mb_per_worker = luigi.IntParameter(default=(16 * 1024))
    split_intervals = luigi.BoolParameter(default=False)
    skip_cleaning = luigi.BoolParameter(default=False)
    log_level = luigi.Parameter(default='WARNING')

    def requires(self):
        logger = logging.getLogger(__name__)
        config = read_config_yml(config_yml_path=self.config_yml_path)
        caller_modes = list()
        if 'callers' in config:
            for k, v in config['callers'].items():
                for c in v:
                    caller_modes.append(f'{k}.{c}')
        annotators = (
            {k for k, v in config['annotators'].items() if v}
            if 'annotators' in config else {'funcotator'}
        )
        common_config = {
            'ref_version': (config.get('reference_version') or 'hg38'),
            'exome': bool(config.get('exome')),
            'memory_mb_per_worker': self.memory_mb_per_worker,
            'n_cpu_per_worker': self.n_cpu_per_worker,
            'gatk_java_options': ' '.join([
                '-Dsamjdk.compression_level=5',
                '-Dsamjdk.use_async_io_read_samtools=true',
                '-Dsamjdk.use_async_io_write_samtools=true',
                '-Dsamjdk.use_async_io_write_tribble=false',
                f'-Xmx{self.memory_mb_per_worker:d}m',
                '-XX:+UseParallelGC',
                f'-XX:ParallelGCThreads={self.n_cpu_per_worker}'
            ]),
            'samtools_memory_per_thread':
            '{:d}M'.format(
                int(self.memory_mb_per_worker / self.n_cpu_per_worker / 20)
            ),
            'save_memory': (self.memory_mb_per_worker < 8 * 1024),
            'split_intervals': self.split_intervals,
            'remove_if_failed': (not self.skip_cleaning),
            **{
                c: fetch_executable(c) for c in {
                    'bcftools', 'bedtools', 'bgzip', 'bwa', 'cutadapt',
                    'fastqc', 'gawk', 'gatk', 'pbzip2', 'pigz', 'samtools',
                    'tabix', 'trim_galore',
                    *(
                        {'python2', 'configureStrelkaSomaticWorkflow.py'}
                        if 'somatic_short_variant.strelka' in caller_modes
                        else set()
                    ),
                    *(
                        {'python3'}
                        if 'germline_short_variant.gatk' in caller_modes
                        else set()
                    ),
                    *(
                        {'python2', 'configureStrelkaGermlineWorkflow.py'}
                        if 'germline_short_variant.strelka' in caller_modes
                        else set()
                    ),
                    *(
                        {'python2', 'configManta.py'}
                        if 'somatic_structual_variant.manta' in caller_modes
                        else set()
                    ),
                    *(
                        {'delly'}
                        if 'somatic_structual_variant.delly' in caller_modes
                        else set()
                    ),
                    *(
                        {
                            'python', 'lumpy', 'lumpyexpress', 'sambamba',
                            'samblaster'
                        } if 'somatic_structual_variant.lumpy' in caller_modes
                        else set()
                    )
                }
            },
            **{
                (k.replace('/', '_') + '_dir_path'): str(
                    Path(self.dest_dir_path).joinpath(k)
                ) for k in {
                    'trim', 'align', 'annotate/bcftools',
                    'annotate/funcotator', 'somatic_snv_indel/gatk',
                    'somatic_snv_indel/strelka', 'germline_snv_indel/gatk',
                    'germline_snv_indel/strelka', 'somatic_sv/manta',
                    'somatic_sv/delly', 'somatic_sv/lumpy', 'somatic_cnv/gatk'
                }
            },
            'ref_dir_path': str(Path(self.ref_dir_path).resolve()),
            'log_dir_path': str(Path(self.log_dir_path).resolve())
        }
        matched_keys = ['tumor', 'normal']
        reference_file_paths = dict([
            ((f'{k}s', [v]) if k == 'ref_fa_path' else (k, v))
            for k, v in self._resolve_input_file_paths(
                path_dict=config['references']
            ).items()
        ])
        task_kwargs = [
            {
                'fq_list': [
                    list(self._resolve_input_file_paths(path_list=r[k]['fq']))
                    for k in matched_keys if r[k].get('fq')
                ],
                'read_groups':
                [(r[k].get('read_group') or dict()) for k in matched_keys],
                'sample_names': [
                    (
                        (r[k].get('read_group') or dict()).get('SM')
                        or parse_fq_id(fq_path=r[k]['fq'][0])
                    ) for k in matched_keys
                ],
                'priority': p, 'cf': common_config,
                'caller_modes': caller_modes, 'annotators': annotators,
                **reference_file_paths
            } for p, r in zip(
                [i * 1000 for i in range(1, (len(config['runs']) + 1))[::-1]],
                config['runs']
            )
        ]
        logger.debug('task_kwargs:' + os.linesep + pformat(task_kwargs))
        return [CallVariants(**d) for d in task_kwargs]

    @staticmethod
    def _resolve_input_file_paths(path_list=None, path_dict=None):
        assert bool(path_list or path_dict)
        if path_list:
            new_list = list()
            for s in path_list:
                p = Path(s).resolve()
                assert p.is_file(), f'file not found: {p}'
                new_list.append(str(p))
            return new_list
        elif path_dict:
            new_dict = dict()
            for k, v in path_dict.items():
                if isinstance(v, str):
                    p = Path(v).resolve()
                    assert p.is_file(), f'file not found: {p}'
                    new_dict[f'{k}_path'] = str(p)
                else:
                    new_dict[f'{k}_paths'] = list()
                    for s in v:
                        p = Path(s).resolve()
                        assert p.is_file(), f'file not found: {p}'
                        new_dict[f'{k}_paths'].append(str(p))
            return new_dict

    def output(self):
        return self.input()

    def run(self):
        logger = logging.getLogger(__name__)
        logger.debug('Task tree:' + os.linesep + deps_tree.print_tree(self))


if __name__ == '__main__':
    luigi.run()
