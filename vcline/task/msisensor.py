#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import create_matched_id
from .align import PrepareBAMNormal, PrepareBAMTumor
from .base import ShellTask
from .ref import (CreateEvaluationIntervalListBED, FetchReferenceFASTA,
                  UncompressBgzipFiles)


@requires(FetchReferenceFASTA)
class ScanMicrosatellites(ShellTask):
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath('msisensor').joinpath(
                Path(self.input()[0].path).stem + '.microsatellites.list'
            )
        )

    def run(self):
        fa_path = self.input()[0].path
        run_id = Path(fa_path).stem
        self.print_log(f'Scan microsatellites:\t{run_id}')
        msisensor = self.cf['msisensor']
        microsatellites_list_path = self.output().path
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=msisensor,
            cwd=Path(microsatellites_list_path).parent,
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {msisensor} scan'
                + f' -d {fa_path} -o {microsatellites_list_path}'
            ),
            input_files_or_dirs=fa_path,
            output_files_or_dirs=microsatellites_list_path
        )


@requires(CreateEvaluationIntervalListBED)
class UncompressEvaluationIntervalListBED(ShellTask):
    cf = luigi.DictParameter()
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.cf['ref_dir_path']).joinpath('msisensor').joinpath(
                Path(self.input()[0].path).stem
            )
        )

    def run(self):
        yield UncompressBgzipFiles(
            bgz_paths=[self.input()[0].path],
            dest_dir_path=str(Path(self.output().path).parent), cf=self.cf
        )


@requires(PrepareBAMTumor, PrepareBAMNormal, ScanMicrosatellites,
          UncompressEvaluationIntervalListBED)
class ScoreMSIWithMSIsensor(ShellTask):
    cf = luigi.DictParameter()
    priority = 20

    def output(self):
        return [
            luigi.LocalTarget(
                Path(self.cf['somatic_msi_msisensor_dir_path']).joinpath(
                    create_matched_id(*[i[0].path for i in self.input()[0:2]])
                    + s
                )
            ) for s in ['', '_dis', '_germline', '_somatic']
        ]

    def run(self):
        output_file_paths = [o.path for o in self.output()]
        run_id = Path(output_file_paths[0]).name
        self.print_log(f'Score MSI with MSIsensor:\t{run_id}')
        msisensor = self.cf['msisensor']
        tumor_bam_path = self.input()[0][0].path
        normal_bam_path = self.input()[1][0].path
        microsatellites_list_path = self.input()[2].path
        bed_path = self.input()[3].path
        output_prefix = Path(output_file_paths[0]).name
        self.setup_shell(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            commands=msisensor, cwd=self.cf['somatic_msi_msisensor_dir_path'],
            remove_if_failed=self.cf['remove_if_failed'],
            quiet=self.cf['quiet']
        )
        self.run_shell(
            args=(
                f'set -e && {msisensor} msi'
                + f' -d {microsatellites_list_path}'
                + f' -n {normal_bam_path}'
                + f' -t {tumor_bam_path}'
                + f' -e {bed_path}'
                + f' -o {output_prefix}'
            ),
            input_files_or_dirs=[
                tumor_bam_path, normal_bam_path, microsatellites_list_path,
                bed_path
            ],
            output_files_or_dirs=output_file_paths
        )


if __name__ == '__main__':
    luigi.run()
