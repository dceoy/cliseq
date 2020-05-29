#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareCRAMNormal, PrepareCRAMTumor
from .base import ShellTask
from .ref import (FetchReferenceFASTA, ScanMicrosatellites,
                  UncompressEvaluationIntervalListBED)
from .samtools import SamtoolsView


@requires(PrepareCRAMTumor, PrepareCRAMNormal, FetchReferenceFASTA,
          ScanMicrosatellites, UncompressEvaluationIntervalListBED)
class ScoreMSIWithMSIsensor(ShellTask):
    cf = luigi.DictParameter()
    priority = 20

    def output(self):
        run_dir = Path(self.cf['somatic_msi_msisensor_dir_path']).joinpath(
            create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(f'{run_dir.name}.msisensor.tsv{s}')
            ) for s in ['', '_dis', '_germline', '_somatic']
        ]

    def run(self):
        output_file_paths = [o.path for o in self.output()]
        run_dir = Path(output_file_paths[0]).parent
        input_targets = yield [
            SamtoolsView(
                input_sam_path=i[0].path,
                output_sam_path=str(
                    run_dir.joinpath(Path(i[0].path).stem + '.bam')
                ),
                fa_path=self.input()[2][0].path,
                samtools=self.cf['samtools'],
                n_cpu=self.cf['n_cpu_per_worker'], remove_input=False,
                index_sam=True, log_dir_path=self.cf['log_dir_path'],
                remove_if_failed=self.cf['remove_if_failed'],
                quiet=self.cf['quiet']
            ) for i in self.input()[0:2]
        ]
        run_id = run_dir.name
        self.print_log(f'Score MSI with MSIsensor:\t{run_id}')
        msisensor = self.cf['msisensor']
        bam_paths = [i[0].path for i in input_targets]
        microsatellites_list_path = self.input()[3].path
        bed_path = self.input()[4].path
        output_path_prefix = Path(output_file_paths[0]).name
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=msisensor, cwd=run_dir,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {msisensor} msi'
                + f' -t {bam_paths[0]} -n {bam_paths[1]}'
                + f' -d {microsatellites_list_path} -e {bed_path}'
                + f' -o {output_path_prefix}'
            ),
            input_files_or_dirs=[
                *bam_paths, microsatellites_list_path, bed_path
            ],
            output_files_or_dirs=[*output_file_paths, run_dir]
        )
        self.run_shell(
            args=('rm -f' + ''.join([f' {p} {p}.bai' for p in bam_paths])),
            input_files_or_dirs=bam_paths
        )


if __name__ == '__main__':
    luigi.run()
