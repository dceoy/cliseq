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
from .funcotator import AnnotateVariantsWithFuncotator


class CallVariants(luigi.WrapperTask):
    ref_fa_paths = luigi.ListParameter()
    fq_list = luigi.ListParameter()
    read_groups = luigi.ListParameter()
    sample_names = luigi.ListParameter()
    dbsnp_vcf_path = luigi.Parameter()
    known_indel_vcf_paths = luigi.ListParameter()
    hapmap_vcf_path = luigi.Parameter()
    gnomad_vcf_path = luigi.Parameter()
    evaluation_interval_path = luigi.Parameter()
    funcotator_somatic_tar_path = luigi.Parameter()
    funcotator_germline_tar_path = luigi.Parameter()
    cf = luigi.DictParameter()
    variant_callers = luigi.DictParameter()
    priority = 100

    def requires(self):
        if any(self.variant_callers.values()):
            return [
                AnnotateVariantsWithFuncotator(
                    variant_caller=k, funcotator_data_source_tar_path=v,
                    ref_fa_paths=self.ref_fa_paths, fq_list=self.fq_list,
                    read_groups=self.read_groups,
                    sample_names=self.sample_names,
                    dbsnp_vcf_path=self.dbsnp_vcf_path,
                    known_indel_vcf_paths=self.known_indel_vcf_paths,
                    gnomad_vcf_path=self.gnomad_vcf_path,
                    hapmap_vcf_path=self.hapmap_vcf_path,
                    evaluation_interval_path=self.evaluation_interval_path,
                    cf=self.cf
                ) for k, v in {
                    'haplotypecaller': self.funcotator_germline_tar_path,
                    'mutect2': self.funcotator_somatic_tar_path,
                    'manta_somatic': self.funcotator_somatic_tar_path,
                    'manta_diploid': self.funcotator_germline_tar_path,
                    'strelka': self.funcotator_somatic_tar_path
                }.items() if self.variant_callers.get(k)
            ]
        else:
            return PrepareCRAMs(
                ref_fa_paths=self.ref_fa_paths, fq_list=self.fq_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                known_indel_vcf_paths=self.known_indel_vcf_paths, cf=self.cf
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
    log_level = luigi.Parameter(default='WARNING')

    def requires(self):
        logger = logging.getLogger(__name__)
        config = read_config_yml(config_yml_path=self.config_yml_path)
        common_config = {
            'ref_version': (config.get('reference_version') or 'hg38'),
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
            **{
                c: fetch_executable(c) for c in [
                    'bgzip', 'bwa', 'cutadapt', 'fastqc', 'gatk', 'pbzip2',
                    'pigz', 'configManta.py',
                    'configureStrelkaSomaticWorkflow.py', 'samtools', 'tabix',
                    'trim_galore'
                ]
            },
            **{
                f'{k}_dir_path': str(Path(self.dest_dir_path).joinpath(k))
                for k in [
                    'trim', 'align', 'haplotypecaller', 'mutect2', 'strelka',
                    'manta'
                ]
            },
            'ref_dir_path': str(Path(self.ref_dir_path).resolve()),
            'log_dir_path': str(Path(self.log_dir_path).resolve())
        }
        fb_keys = ['foreground', 'background']
        reference_file_paths = self._resolve_input_file_paths(
            path_dict=config['references']
        )
        variant_callers = (
            config.get('variant_callers') or {
                'haplotypecaller': True, 'mutect2': True,
                'strelka': True, 'manta': True
            }
        )
        task_kwargs = [
            {
                'fq_list': [
                    list(self._resolve_input_file_paths(path_list=r[k]['fq']))
                    for k in fb_keys if r[k].get('fq')
                ],
                'read_groups':
                [(r[k].get('read_group') or dict()) for k in fb_keys],
                'sample_names': [
                    (
                        (r[k].get('read_group') or dict()).get('SM')
                        or parse_fq_id(fq_path=r[k]['fq'][0])
                    ) for k in fb_keys
                ],
                'cf': common_config, 'variant_callers': variant_callers,
                **reference_file_paths
             } for r in config['runs']
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
