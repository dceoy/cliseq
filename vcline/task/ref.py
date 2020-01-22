#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .base import ShellTask


class FetchGenomeFASTA(ShellTask):
    ref_fa_list = luigi.ListParameter()
    cf = luigi.DictParameter()
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            self.ref_fa_list[0]['src'] if (
                len(self.ref_fa_list) == 1
                and (not self.ref_fa_list[0]['is_url'])
                and (str(Path(self.ref_fa_list[0]['src']).parent)
                     == self.cf['ref_dir_path'])
                and self.ref_fa_list[0]['src'].endswith(('.fa', '.fasta'))
            ) else str(
                Path(self.cf['ref_dir_path']).joinpath(
                    '.'.join([
                        Path(Path(d['src']).stem).stem
                        for d in self.ref_fa_list
                    ]) + '.fa'
                )
            )
        )

    def run(self):
        fa_path = self.output().path
        run_id = Path(fa_path).stem
        print_log(f'Create a reference FASTA:\t{run_id}')
        cat = self.cf['cat']
        curl = self.cf['curl']
        pigz = self.cf['pigz']
        pbzip2 = self.cf['pbzip2']
        n_cpu = self.cf['n_cpu_per_worker']
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        args = [
            f'{cat} --version',
            f'{curl} --version',
            f'{pigz} --version',
            f'{pbzip2} --version'
        ]
        input_files = list()
        for i, d in enumerate(self.ref_fa_list):
            s = d['src']
            if not d['is_url']:
                input_files.append(s)
                if s.endswith('.gz'):
                    a = f'{cat} {s} | {pigz} -p {n_cpu} -dc -'
                elif s.endswith('.bz2'):
                    a = f'{cat} {s} | {pbzip2} -p# {n_cpu} -dc -'
                else:
                    a = '{cat} {s}'
            elif s.endswith('.gz'):
                a = f'{curl} -LS {s} | {pigz} -p {n_cpu} -dc -'
            elif s.endswith('.bz2'):
                a = f'{curl} -LS {s} | {pbzip2} -p# {n_cpu} -dc -'
            else:
                a = f'{curl} -LS {s}'
            r = '>' if i == 0 else '>>'
            args.append(f'set -eo pipefail && {a} {r} {fa_path}')
        self.run_bash(args=args, input_files=input_files, output_files=fa_path)


@requires(FetchGenomeFASTA)
class CreateFASTAIndex(ShellTask):
    cf = luigi.DictParameter()
    priority = 8

    def output(self):
        return luigi.LocalTarget(f'{self.input().path}.fai')

    def run(self):
        fa_path = self.input().path
        run_id = Path(fa_path).stem
        print_log(f'Create a FASTA index:\t{run_id}')
        samtools = self.cf['samtools']
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        self.run_bash(
            args=[
                f'{samtools} 2>&1 | grep -e "Version:"',
                f'set -e && {samtools} faidx {fa_path}'
            ],
            input_files=fa_path, output_files=self.output().path
        )


@requires(FetchGenomeFASTA)
class CreateBWAIndexes(ShellTask):
    cf = luigi.DictParameter()
    priority = 9

    def output(self):
        return [
            luigi.LocalTarget(f'{self.input().path}.{s}')
            for s in ['pac', 'bwt', 'ann', 'amb', 'sa']
        ]

    def run(self):
        fa_path = self.input().path
        run_id = Path(fa_path).stem
        print_log(f'Create BWA indexes:\t{run_id}')
        bwa = self.cf['bwa']
        self.setup_bash(
            run_id=run_id, log_dir_path=self.cf['log_dir_path'],
            work_dir_path=self.cf['ref_dir_path']
        )
        self.run_bash(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'set -e && {bwa} index {fa_path}'
            ],
            input_files=fa_path,
            output_files=[o.path for o in self.output()]
        )


if __name__ == '__main__':
    luigi.run()
