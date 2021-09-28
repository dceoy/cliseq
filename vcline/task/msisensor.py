#!/usr/bin/env python

from pathlib import Path

import luigi
from ftarc.task.resource import FetchReferenceFasta
from ftarc.task.samtools import SamtoolsView
from luigi.util import requires

from .core import VclineTask
from .cram import PrepareCramNormal, PrepareCramTumor
from .resource import CreateEvaluationIntervalListBed


@requires(FetchReferenceFasta)
class ScanMicrosatellites(VclineTask):
    cf = luigi.DictParameter()
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        fa = Path(self.input()[0].path).resolve()
        return luigi.LocalTarget(
            fa.parent.joinpath(f'{fa.stem}.microsatellites.tsv')
        )

    def run(self):
        fa = Path(self.input()[0].path)
        run_id = fa.stem
        self.print_log(f'Scan microsatellites:\t{run_id}')
        msisensor = self.cf['msisensor']
        output_tsv = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=msisensor, cwd=fa.parent,
            **self.sh_config
        )
        self.run_shell(
            args=f'set -e && {msisensor} scan -d {fa} -o {output_tsv}',
            input_files_or_dirs=fa, output_files_or_dirs=output_tsv
        )


@requires(CreateEvaluationIntervalListBed)
class UncompressEvaluationIntervalListBed(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        bed_gz = Path(self.input()[0].path)
        return luigi.LocalTarget(bed_gz.parent.joinpath(bed_gz.stem))

    def run(self):
        output_bed = Path(self.output().path)
        run_id = output_bed.stem
        self.print_log(f'Uncompress bgzip files:\t{run_id}')
        bed_gz = Path(self.input()[0].path)
        bgzip = self.cf['bgzip']
        self.setup_shell(
            run_id=run_id, commands=bgzip, cwd=output_bed.parent,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {bgzip} -@ {self.n_cpu} -dc {bed_gz}'
                + f' > {output_bed}'
            ),
            input_files_or_dirs=bed_gz, output_files_or_dirs=output_bed
        )


@requires(PrepareCramTumor, PrepareCramNormal, FetchReferenceFasta,
          ScanMicrosatellites, UncompressEvaluationIntervalListBed)
class ScoreMsiWithMsisensor(VclineTask):
    cf = luigi.DictParameter()
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 20

    def output(self):
        run_dir = Path(self.cf['somatic_msi_msisensor_dir_path']).joinpath(
            self.create_matched_id(*[i[0].path for i in self.input()[0:2]])
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(f'{run_dir.name}.msisensor.tsv{s}')
            ) for s in ['', '_dis', '_germline', '_somatic']
        ]

    def run(self):
        output_files = [Path(o.path) for o in self.output()]
        run_dir = output_files[0].parent
        input_targets = yield [
            SamtoolsView(
                input_sam_path=i[0].path,
                output_sam_path=str(
                    run_dir.joinpath(Path(i[0].path).stem + '.bam')
                ),
                fa_path=self.input()[2][0].path,
                samtools=self.cf['samtools'], n_cpu=self.n_cpu,
                remove_input=False, index_sam=True, sh_config=self.sh_config
            ) for i in self.input()[0:2]
        ]
        run_id = run_dir.name
        self.print_log(f'Score MSI with MSIsensor:\t{run_id}')
        msisensor = self.cf['msisensor']
        bams = [Path(i[0].path) for i in input_targets]
        microsatellites_list = Path(self.input()[3].path)
        bed = Path(self.input()[4].path)
        output_path_prefix = output_files[0].name
        self.setup_shell(
            run_id=run_id, commands=msisensor, cwd=run_dir,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {msisensor} msi'
                + f' -t {bams[0]} -n {bams[1]}'
                + f' -d {microsatellites_list} -e {bed}'
                + f' -o {output_path_prefix}'
            ),
            input_files_or_dirs=[*bams, microsatellites_list, bed],
            output_files_or_dirs=[*output_files, run_dir]
        )
        self.remove_files_and_dirs(*bams, *[f'{p}.bai' for p in bams])


if __name__ == '__main__':
    luigi.run()
