#!/usr/bin/env python

from pathlib import Path

import luigi

from .core import VclineTask


class CallSomaticCnvWithCnvkit(VclineTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    refflat_txt_path = luigi.Parameter()
    access_bed_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    cnvkitpy = luigi.Parameter(default='cnvkit.py')
    samtools = luigi.Parameter(default='samtools')
    rscript = luigi.Parameter(default='Rscript')
    seq_method = luigi.Parameter(default='wgs')
    n_cpu = luigi.IntParameter(default=1)
    drop_low_coverage = luigi.BoolParameter(default=True)
    short_names = luigi.BoolParameter(default=True)
    diagram = luigi.BoolParameter(default=True)
    scatter = luigi.BoolParameter(default=True)
    sh_config = luigi.DictParameter(default=dict())
    priority = 20

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        tumor_stem = Path(self.tumor_cram_path).stem
        normal_stem = Path(self.normal_cram_path).stem
        access_stem = Path(self.access_bed_path).stem
        return [
            luigi.LocalTarget(run_dir.joinpath(n)) for n in (
                [
                    (tumor_stem + s) for s in (
                        [
                            '.call.seg', '.seg', '.call.cns', '.cns',
                            '.bintest.cns', '.cnr', '.targetcoverage.cnn',
                            '.antitargetcoverage.cnn'
                        ] + (['.diagram.pdf'] if self.diagram else list())
                        + (['.scatter.pdf'] if self.scatter else list())
                    )
                ] + [
                    (normal_stem + s) for s in [
                        '.targetcoverage.cnn', '.antitargetcoverage.cnn',
                        '.reference.cnn'
                    ]
                ] + [
                    (access_stem + s) for s in [
                        '.target.bed', '.antitarget.bed'
                    ]
                ]
            )
        ]

    def run(self):
        run_id = self.create_matched_id(
            self.tumor_cram_path, self.normal_cram_path
        )
        self.print_log(f'Call somatic CNVs with CNVkit:\t{run_id}')
        tumor_cram = Path(self.tumor_cram_path).resolve()
        normal_cram = Path(self.normal_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        access_bed = Path(self.access_bed_path).resolve()
        refflat_txt = Path(self.refflat_txt_path).resolve()
        output_files = [Path(o.path) for o in self.output()]
        run_dir = output_files[0].parent
        output_ref_cnn = run_dir.joinpath(f'{normal_cram.stem}.reference.cnn')
        output_call_cns = output_files[2]
        output_cns = output_files[3]
        self.setup_shell(
            run_id=run_id,
            commands=[self.cnvkitpy, self.samtools, self.rscript], cwd=run_dir,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.cnvkitpy} batch'
                + f' --seq-method={self.seq_method}'
                + f' --fasta={fa}'
                + f' --access={access_bed}'
                + f' --annotate={refflat_txt}'
                + f' --processes={self.n_cpu}'
                + (' --drop-low-coverage' if self.drop_low_coverage else '')
                + (' --short-names' if self.short_names else '')
                + f' --output-dir={run_dir}'
                + f' --output-reference={output_ref_cnn}'
                + f' --normal={normal_cram}'
                + f' {tumor_cram}'
            ),
            input_files_or_dirs=[
                tumor_cram, normal_cram, fa, access_bed, refflat_txt
            ],
            output_files_or_dirs=[
                *[f for f in output_files[2:] if f.suffix != '.pdf'], run_dir
            ]
        )
        for o in [output_call_cns, output_cns]:
            output_seg = run_dir.joinpath(f'{o.stem}.seg')
            self.run_shell(
                args=(
                    f'set -e && {self.cnvkitpy} export seg'
                    + f' --output={output_seg} {o}'
                ),
                input_files_or_dirs=o, output_files_or_dirs=output_seg
            )
        for c in ['diagram', 'scatter']:
            if getattr(self, c):
                graph_pdf = run_dir.joinpath(f'{output_cns.stem}.{c}.pdf')
                self.run_shell(
                    args=(
                        f'set -e && {self.cnvkitpy} {c}'
                        + f' --output={graph_pdf} {output_cns}'
                    ),
                    input_files_or_dirs=output_cns,
                    output_files_or_dirs=graph_pdf
                )


if __name__ == '__main__':
    luigi.run()
