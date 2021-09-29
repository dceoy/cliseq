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
    sh_config = luigi.DictParameter(default=dict())
    priority = 20

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        tumor_stem = Path(self.tumor_cram_path).stem
        normal_stem = Path(self.normal_cram_path).stem
        return [
            luigi.LocalTarget(run_dir.joinpath(n)) for n in (
                [
                    (tumor_stem + s) for s in [
                        '.call.seg', '.call.cns', '.cns', '.bintest.cns',
                        '.cnr', '.targetcoverage.cnn',
                        '.antitargetcoverage.cnn', '-diagram.pdf',
                        '-scatter.png'
                    ]
                ] + [
                    (normal_stem + s) for s in [
                        '.targetcoverage.cnn', '.antitargetcoverage.cnn',
                        '.reference.cnn'
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
        output_call_cns = output_files[0]
        run_dir = output_call_cns.parent
        output_ref_cnn = run_dir.joinpath(f'{normal_cram.stem}.reference.cnn')
        output_call_seg = run_dir.joinpath(f'{output_call_cns.stem}.seg')
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
                + ' --drop-low-coverage --diagram --scatter'
                + f' --output-dir={run_dir}'
                + f' --output-reference={output_ref_cnn}'
                + f' --normal={normal_cram}'
                + f' {tumor_cram}'
            ),
            input_files_or_dirs=[
                tumor_cram, normal_cram, fa, access_bed, refflat_txt
            ],
            output_files_or_dirs=output_files[1:]
        )
        self.run_shell(
            args=(
                f'set -e && {self.cnvkitpy} export seg'
                + f' --output={output_call_seg}'
                + f' {output_call_cns}'
            ),
            input_files_or_dirs=output_call_cns,
            output_files_or_dirs=output_call_seg
        )


if __name__ == '__main__':
    luigi.run()
