#!/usr/bin/env python

import logging
import os
import re
import sys
from itertools import chain
from pathlib import Path
from pprint import pformat

import luigi
import yaml
from luigi.tools import deps_tree

from ..cli.util import (create_matched_id, fetch_executable, parse_cram_id,
                        parse_fq_id, read_yml)
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .base import BaseTask, ShellTask
from .callcopyratiosegments import CallCopyRatioSegmentsMatched
from .canvas import CallSomaticCopyNumberVariantsWithCanvas
from .delly import CallStructualVariantsWithDelly
from .funcotator import FuncotateSegments, FuncotateVariants
from .haplotypecaller import FilterVariantTranches
from .lumpy import CallStructualVariantsWithLumpy
from .manta import CallStructualVariantsWithManta
from .msisensor import ScoreMSIWithMSIsensor
from .mutect2 import FilterMutectCalls
from .snpeff import AnnotateVariantsWithSnpEff
from .strelka import (CallGermlineVariantsWithStrelka,
                      CallSomaticVariantsWithStrelka)


class RunVariantCaller(BaseTask):
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
    funcotator_data_src_tar_path = luigi.Parameter()
    snpeff_config_path = luigi.Parameter()
    cf = luigi.DictParameter()
    caller_mode = luigi.ListParameter()
    annotators = luigi.ListParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = 10

    def requires(self):
        if 'germline_short_variant.gatk' == self.caller_mode:
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
        elif 'somatic_short_variant.gatk' == self.caller_mode:
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
        elif 'somatic_structual_variant.manta' == self.caller_mode:
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
        elif 'somatic_short_variant.strelka' == self.caller_mode:
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
        elif 'germline_short_variant.strelka' == self.caller_mode:
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
        elif 'somatic_structual_variant.delly' == self.caller_mode:
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
        elif 'somatic_structual_variant.lumpy' == self.caller_mode:
            return CallStructualVariantsWithLumpy(
                fq_list=self.fq_list, cram_list=self.cram_list,
                read_groups=self.read_groups, sample_names=self.sample_names,
                ref_fa_path=self.ref_fa_path,
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                mills_indel_vcf_path=self.mills_indel_vcf_path,
                known_indel_vcf_path=self.known_indel_vcf_path,
                evaluation_interval_path=self.evaluation_interval_path,
                cf=self.cf
            )
        elif 'somatic_copy_number_variation.gatk' == self.caller_mode:
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
        elif 'somatic_copy_number_variation.canvas' == self.caller_mode:
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
        elif 'somatic_msi.msisensor' == self.caller_mode:
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
            raise ValueError(f'invalid caller mode: {self.caller_mode}')

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
        if 'somatic_structual_variant.delly' == self.caller_mode:
            suffix_dict = {'funcotator': None, 'snpeff': '.vcf.gz'}
        elif 'somatic_structual_variant.manta' == self.caller_mode:
            suffix_dict = {
                'funcotator': '.manta.somaticSV.vcf.gz',
                'snpeff': ('.manta.somaticSV.vcf.gz', '.diploidSV.vcf.gz')
            }
        elif 'germline_short_variant.strelka' == self.caller_mode:
            suffix_dict = {
                k: '.strelka.germline.variants.vcf.gz'
                for k in ['funcotator', 'snpeff']
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
        for k, v in self._find_annotation_targets().items():
            if k == 'funcotator':
                common_kwargs = {
                    'data_src_tar_path': self.funcotator_data_src_tar_path,
                    'ref_fa_path': self.ref_fa_path, 'cf': self.cf
                }
                yield [
                    (
                        FuncotateSegments(input_seg_path=p, **common_kwargs)
                        if p.endswith('.seg')
                        else FuncotateVariants(
                            input_vcf_path=p, normalize_vcf=self.normalize_vcf,
                            **common_kwargs
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


class CallVariants(ShellTask):
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
    caller_modes = luigi.ListParameter()
    annotators = luigi.ListParameter()
    normalize_vcf = luigi.BoolParameter(default=True)
    priority = luigi.IntParameter(default=100)

    def requires(self):
        common_kwargs = {
            'fq_list': self.fq_list, 'read_groups': self.read_groups,
            'cram_list': self.cram_list, 'sample_names': self.sample_names,
            'ref_fa_path': self.ref_fa_path,
            'dbsnp_vcf_path': self.dbsnp_vcf_path,
            'mills_indel_vcf_path': self.mills_indel_vcf_path,
            'known_indel_vcf_path': self.known_indel_vcf_path, 'cf': self.cf
        }
        return [
            PrepareCRAMTumor(**common_kwargs),
            PrepareCRAMNormal(**common_kwargs),
            *[
                RunVariantCaller(
                    hapmap_vcf_path=self.hapmap_vcf_path,
                    gnomad_vcf_path=self.gnomad_vcf_path,
                    evaluation_interval_path=self.evaluation_interval_path,
                    cnv_black_list_path=self.cnv_black_list_path,
                    genomesize_xml_path=self.genomesize_xml_path,
                    kmer_fa_path=self.kmer_fa_path,
                    exome_manifest_path=self.exome_manifest_path,
                    funcotator_data_src_tar_path=(
                        self.funcotator_germline_tar_path
                        if m.startswith('.germline') else
                        self.funcotator_somatic_tar_path
                    ),
                    snpeff_config_path=self.snpeff_config_path, caller_mode=m,
                    annotators=self.annotators,
                    normalize_vcf=self.normalize_vcf, **common_kwargs
                ) for m in self.caller_modes
            ]
        ]

    def output(self):
        return self.input()

    def run(self):
        bam_paths = [
            str(
                Path(self.cf['align_dir_path']).joinpath(
                    Path(i[0].path).stem + '.bam'
                )
            ) for i in self.input()[:2]
        ]
        if any([Path(p).is_file() for p in bam_paths]):
            run_id = create_matched_id(*bam_paths)
            self.print_log(f'Remove BAM files:\t{run_id}')
            self.setup_shell(
                run_id=run_id, log_dir_path=self.cf['log_dir_path'],
                cwd=self.cf['align_dir_path'],
                remove_if_failed=self.cf['remove_if_failed']
            )
            self.run_shell(
                args=('rm -f' + ''.join([f' {p}*' for p in bam_paths]))
            )


class PrintEnvVersions(ShellTask):
    cf = luigi.DictParameter()
    priority = sys.maxsize

    def run(self):
        python = sys.executable
        self.print_log(f'Print environment versions: {python}')
        self.setup_shell(
            run_id='env', log_dir_path=self.cf['log_dir_path'],
            commands=[python, self.cf['java']]
        )
        self.run_shell(
            args=[f'{python} -m pip --version', f'{python} -m pip freeze']
        )

    def complete(self):
        return True


class RunAnalyticalPipeline(BaseTask):
    config_yml_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    log_dir_path = luigi.Parameter(default='log')
    ref_dir_path = luigi.Parameter(default='')
    n_cpu_per_worker = luigi.IntParameter(default=1)
    memory_mb_per_worker = luigi.IntParameter(default=(16 * 1024))
    skip_cleaning = luigi.BoolParameter(default=False)
    log_level = luigi.Parameter(default='WARNING')

    def requires(self):
        logger = logging.getLogger(__name__)
        config = self._read_config_yml(config_yml_path=self.config_yml_path)
        caller_modes = list()
        if 'callers' in config:
            for k, v in config['callers'].items():
                for t, b in v.items():
                    if b:
                        caller_modes.append(f'{k}.{t}')
        if ({m for m in caller_modes if m.endswith('.canvas')}
                and self.n_cpu_per_worker >= 8
                and self.memory_mb_per_worker >= 25 * 1024):
            logger.warning('Canvas requires 8 CPUs and 25 GB RAM.')
        annotators = (
            {k for k, v in config['annotators'].items() if v}
            if 'annotators' in config else {'funcotator', 'snpeff'}
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
            'remove_if_failed': (not self.skip_cleaning),
            **{
                c: fetch_executable(c) for c in {
                    'bcftools', 'bedtools', 'bgzip', 'bwa', 'cutadapt',
                    'fastqc', 'gawk', 'gatk', 'java', 'pbzip2', 'pigz',
                    'samtools', 'tabix', 'trim_galore',
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
                    ),
                    *(
                        {'R'} if (
                            'somatic_copy_number_variation.gatk'
                            in caller_modes
                        ) else set()
                    ),
                    *(
                        {'Canvas', 'dotnet'} if (
                            'somatic_copy_number_variation.canvas'
                            in caller_modes
                        ) else set()
                    ),
                    *({'snpEff'} if 'snpeff' in annotators else set()),
                    *(
                        {'msisensor'}
                        if 'somatic_msi.msisensor' in caller_modes else set()
                    )
                }
            },
            **{
                (k.replace('/', '_') + '_dir_path'): str(
                    Path(self.dest_dir_path).joinpath(k)
                ) for k in {
                    'trim', 'align', 'postproc/bcftools',
                    'postproc/funcotator', 'postproc/snpeff',
                    'somatic_snv_indel/gatk', 'somatic_snv_indel/strelka',
                    'germline_snv_indel/gatk', 'germline_snv_indel/strelka',
                    'somatic_sv/manta', 'somatic_sv/delly', 'somatic_sv/lumpy',
                    'somatic_cnv/gatk', 'somatic_cnv/canvas',
                    'somatic_msi/msisensor'
                }
            },
            'ref_dir_path': str(Path(self.ref_dir_path).resolve()),
            'log_dir_path': str(Path(self.log_dir_path).resolve())
        }
        reference_file_paths = self._resolve_input_file_paths(
            path_dict=config['references']
        )
        task_kwargs = [
            {
                'priority': p, 'cf': common_config,
                'caller_modes': caller_modes, 'annotators': annotators,
                **self._determine_input_samples(run_dict=r),
                **reference_file_paths
            } for p, r in zip(
                [i * 1000 for i in range(1, (len(config['runs']) + 1))[::-1]],
                config['runs']
            )
        ]
        yaml.dump({
            'SAMPLES': [
                dict(zip(['TUMOR', 'NORMAL'], d['sample_names']))
                for d in task_kwargs
            ]
        })
        logger.info('task_kwargs:' + os.linesep + pformat(task_kwargs))
        return [
            PrintEnvVersions(cf=common_config),
            *[CallVariants(**d) for d in task_kwargs]
        ]

    def _read_config_yml(self, config_yml_path):
        config = read_yml(path=str(Path(config_yml_path).resolve()))
        assert isinstance(config, dict), config
        for k in ['references', 'runs']:
            assert config.get(k), k
        assert isinstance(config['references'], dict), config['references']
        for k in ['ref_fa', 'dbsnp_vcf', 'mills_indel_vcf', 'known_indel_vcf',
                  'hapmap_vcf', 'gnomad_vcf', 'evaluation_interval',
                  'funcotator_germline_tar', 'funcotator_somatic_tar',
                  'snpeff_config', 'exome_manifest']:
            v = config['references'].get(k)
            if k == 'ref_fa' and isinstance(v, list) and v:
                assert self._has_unique_elements(v), k
                for s in v:
                    assert isinstance(s, str), k
            elif v:
                assert isinstance(v, str), k
        assert isinstance(config['runs'], list), config['runs']
        for r in config['runs']:
            assert isinstance(r, dict), r
            assert set(r.keys()).intersection({'tumor', 'normal'}), r
            for t in ['tumor', 'normal']:
                s = r[t]
                assert isinstance(s, dict), s
                assert (s.get('fq') or s.get('cram')), s
                if s.get('fq'):
                    assert isinstance(s['fq'], list), s
                    assert self._has_unique_elements(s['fq']), s
                    assert (len(s['fq']) <= 2), s
                else:
                    assert isinstance(s['cram'], str), s
                if s.get('read_group'):
                    assert isinstance(s['read_group'], dict), s
                    for k, v in s['read_group'].items():
                        assert re.fullmatch(r'[A-Z]{2}', k), k
                        assert isinstance(v, str), k
                if s.get('sample_name'):
                    assert isinstance(s['sample_name'], str), s
        return config

    @staticmethod
    def _has_unique_elements(elements):
        return len(set(elements)) == len(tuple(elements))

    @staticmethod
    def _resolve_file_path(path):
        p = Path(path).resolve()
        assert p.is_file(), f'file not found: {p}'
        return str(p)

    def _resolve_input_file_paths(self, path_list=None, path_dict=None):
        if path_list:
            return [self._resolve_file_path(s) for s in path_list]
        elif path_dict:
            new_dict = dict()
            for k, v in path_dict.items():
                if isinstance(v, str):
                    new_dict[f'{k}_path'] = self._resolve_file_path(v)
                elif v:
                    new_dict[f'{k}_paths'] = [
                        self._resolve_file_path(s) for s in v
                    ]
            return new_dict

    def _determine_input_samples(self, run_dict):
        tn = [run_dict[i] for i in ['tumor', 'normal']]
        cram_paths = [d.get('cram') for d in tn]
        if all(cram_paths):
            cram_list = self._resolve_input_file_paths(path_list=cram_paths)
            return {
                'fq_list': list(), 'read_groups': list(),
                'cram_list': cram_list,
                'sample_names': [
                    (
                        d.get('sample_name')
                        or parse_cram_id(cram_path=d['cram'])
                    ) for d in tn
                ]
            }
        else:
            return {
                'fq_list': [
                    list(self._resolve_input_file_paths(path_list=d['fq']))
                    for d in tn
                ],
                'read_groups': [(d.get('read_group') or dict()) for d in tn],
                'cram_list': list(),
                'sample_names': [
                    (
                        (d.get('read_group') or dict()).get('SM')
                        or parse_fq_id(fq_path=d['fq'][0])
                    ) for d in tn
                ]
            }

    def output(self):
        return self.input()

    def run(self):
        logger = logging.getLogger(__name__)
        logger.debug('Task tree:' + os.linesep + deps_tree.print_tree(self))


if __name__ == '__main__':
    luigi.run()
