#!/usr/bin/env python

from itertools import product
from math import floor
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMs
from .base import ShellTask
from .manta import CallStructualVariantsWithManta
from .ref import (CreateEvaluationIntervalListBED, CreateFASTAIndex,
                  FetchReferenceFASTA)


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          CreateEvaluationIntervalListBED, CallStructualVariantsWithManta)
class CallSomaticVariantsWithStrelka(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['strelka_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0]]
                        )
                    ).joinpath(f'results/variants/somatic.{v}{i}')
                )
            )
            for v, i in product(['snvs.vcf.gz', 'indels.vcf.gz'], ['', '.tbi'])
        ]

    def run(self):
        output_file_paths = [o.path for o in self.output()]
        run_dir_path = str(Path(output_file_paths[0]).parent.parent.parent)
        run_id = Path(run_dir_path).name
        self.print_log(f'Call somatic variants with Strelka:\t{run_id}')
        config_script = self.cf['configureStrelkaSomaticWorkflow.py']
        run_script = str(Path(run_dir_path).joinpath('runWorkflow.py'))
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 1)
        input_cram_paths = [i[0].path for i in self.input()[0]]
        fa_path = self.input()[1].path
        fai_path = self.input()[2].path
        bed_path = self.input()[3][0].path
        manta_indel_vcf_path = str(
            Path(self.input()[4][0].path).parent.joinpath(
                'candidateSmallIndels.vcf.gz'
            )
        )
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=config_script, cwd=self.cf['strelka_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {config_script}'
                + f' --tumorBam={input_cram_paths[0]}'
                + f' --normalBam={input_cram_paths[1]}'
                + f' --referenceFasta={fa_path}'
                + f' --indelCandidates={manta_indel_vcf_path}'
                + f' --callRegions={bed_path}'
                + f' --runDir={run_dir_path}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[
                *input_cram_paths, fa_path, fai_path, manta_indel_vcf_path,
                bed_path
            ],
            output_files_or_dirs=[run_script, run_dir_path]
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {run_script}'
                + f' --jobs={n_cpu}'
                + f' --memGb={memory_gb}'
                + ' --mode=local'
            ),
            input_files_or_dirs=[
                run_script, *input_cram_paths, fa_path, fai_path,
                manta_indel_vcf_path, bed_path
            ],
            output_files_or_dirs=[*output_file_paths, run_dir_path]
        )


@requires(PrepareCRAMs, FetchReferenceFASTA, CreateFASTAIndex,
          CreateEvaluationIntervalListBED)
class CallGermlineVariantsWithStrelka(ShellTask):
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['strelka_dir_path']).joinpath(
                        Path(self.input()[0][1][0].path).stem
                    ).joinpath(f'results/variants/{v}{i}')
                )
            ) for v, i
            in product(['genome.vcf.gz', 'variants.vcf.gz'], ['', '.tbi'])
        ]

    def run(self):
        output_file_paths = [o.path for o in self.output()]
        run_dir_path = str(Path(output_file_paths[0]).parent.parent.parent)
        run_id = Path(run_dir_path).name
        self.print_log(f'Call germline variants with Strelka:\t{run_id}')
        config_script = self.cf['configureStrelkaGermlineWorkflow.py']
        run_script = str(Path(run_dir_path).joinpath('runWorkflow.py'))
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 1)
        input_cram_path = self.input()[0][1][0].path
        fa_path = self.input()[1].path
        fai_path = self.input()[2].path
        bed_path = self.input()[3][0].path
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=config_script, cwd=self.cf['strelka_dir_path'],
            remove_if_failed=self.cf['remove_if_failed']
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {config_script}'
                + f' --bam={input_cram_path}'
                + f' --referenceFasta={fa_path}'
                + f' --callRegions={bed_path}'
                + f' --runDir={run_dir_path}'
                + (' --exome' if self.cf['exome'] else '')
            ),
            input_files_or_dirs=[input_cram_path, fa_path, fai_path, bed_path],
            output_files_or_dirs=[run_script, run_dir_path]
        )
        self.run_shell(
            args=(
                f'set -e && PYTHONPATH="{pythonpath}" && {run_script}'
                + f' --jobs={n_cpu}'
                + f' --memGb={memory_gb}'
                + ' --mode=local'
            ),
            input_files_or_dirs=[
                run_script, input_cram_path, fa_path, fai_path, bed_path
            ],
            output_files_or_dirs=[*output_file_paths, run_dir_path]
        )


if __name__ == '__main__':
    luigi.run()
