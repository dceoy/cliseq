#!/usr/bin/env python

from itertools import product
from math import floor
from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareNormalCRAM, PrepareTumorCRAM
from .base import ShellTask
from .manta import CallStructualVariantsWithManta
from .ref import (CreateEvaluationIntervalListBED, CreateFASTAIndex,
                  FetchReferenceFASTA)


@requires(PrepareTumorCRAM, PrepareNormalCRAM, FetchReferenceFASTA,
          CreateFASTAIndex, CreateEvaluationIntervalListBED,
          CallStructualVariantsWithManta)
class CallSomaticVariantsWithStrelka(ShellTask):
    cf = luigi.DictParameter()
    result_file_names = [
        (v + s) for v, s in
        product(['somatic.snvs.vcf.gz', 'somatic.indels.vcf.gz'], ['', '.tbi'])
    ]
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['somatic_sv_strelka_dir_path']).joinpath(
                        create_matched_id(
                            *[i[0].path for i in self.input()[0:2]]
                        ) + f'.strelka.{n}'
                    )
                )
            ) for n in self.result_file_names
        ]

    def run(self):
        output_link_paths = [o.path for o in self.output()]
        run_id = '.'.join(Path(output_link_paths[0]).name.split('.')[:-5])
        self.print_log(f'Call somatic variants with Strelka:\t{run_id}')
        config_script = self.cf['configureStrelkaSomaticWorkflow.py']
        root_dir_path = self.cf['somatic_sv_strelka_dir_path']
        run_dir_path = str(Path(root_dir_path).joinpath(run_id))
        run_script = str(Path(run_dir_path).joinpath('runWorkflow.py'))
        python2 = self.cf['python2']
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 1)
        input_cram_paths = [i[0].path for i in self.input()[0:2]]
        fa_path = self.input()[2].path
        fai_path = self.input()[3].path
        bed_path = self.input()[4][0].path
        manta_indel_vcf_path = [
            i.path for i in self.input()[5]
            if i.path.endswith('candidateSmallIndels.vcf.gz')
        ][0]
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        result_file_paths = [
            str(Path(run_dir_path).joinpath(f'results/variants/{n}'))
            for n in self.result_file_names
        ]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=[python2, config_script], cwd=root_dir_path,
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
            output_files_or_dirs=[*result_file_paths, run_dir_path]
        )
        self.run_shell(
            args=[
                'ln -s {0} {1}'.format(Path(i).relative_to(root_dir_path), o)
                for i, o in zip(result_file_paths, output_link_paths)
            ],
            input_files_or_dirs=result_file_paths,
            output_files_or_dirs=output_link_paths
        )


@requires(PrepareNormalCRAM, FetchReferenceFASTA, CreateFASTAIndex,
          CreateEvaluationIntervalListBED)
class CallGermlineVariantsWithStrelka(ShellTask):
    cf = luigi.DictParameter()
    result_file_names = [
        (v + s) for v, s
        in product(['genome.vcf.gz', 'variants.vcf.gz'], ['', '.tbi'])
    ]
    priority = 10

    def output(self):
        return [
            luigi.LocalTarget(
                str(
                    Path(self.cf['germline_sv_strelka_dir_path']).joinpath(
                        Path(self.input()[0][0].path).stem
                        + f'.strelka.germline.{n}'
                    )
                )
            ) for n in self.result_file_names
        ]

    def run(self):
        output_link_paths = [o.path for o in self.output()]
        run_id = '.'.join(Path(output_link_paths[0]).name.split('.')[:-5])
        self.print_log(f'Call germline variants with Strelka:\t{run_id}')
        config_script = self.cf['configureStrelkaGermlineWorkflow.py']
        root_dir_path = self.cf['germline_sv_strelka_dir_path']
        run_dir_path = str(Path(root_dir_path).joinpath(run_id))
        run_script = str(Path(run_dir_path).joinpath('runWorkflow.py'))
        n_cpu = self.cf['n_cpu_per_worker']
        memory_gb = max(floor(self.cf['memory_mb_per_worker'] / 1024), 1)
        input_cram_path = self.input()[0][0].path
        fa_path = self.input()[1].path
        fai_path = self.input()[2].path
        bed_path = self.input()[3][0].path
        pythonpath = Path(config_script).parent.parent.joinpath('lib/python')
        result_file_paths = [
            str(Path(run_dir_path).joinpath(f'results/variants/{n}'))
            for n in self.result_file_names
        ]
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=config_script, cwd=root_dir_path,
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
            output_files_or_dirs=[*result_file_paths, run_dir_path]
        )
        self.run_shell(
            args=[
                'ln -s {0} {1}'.format(Path(i).relative_to(root_dir_path), o)
                for i, o in zip(result_file_paths, output_link_paths)
            ],
            input_files_or_dirs=result_file_paths,
            output_files_or_dirs=output_link_paths
        )


if __name__ == '__main__':
    luigi.run()
