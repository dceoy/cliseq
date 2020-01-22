#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import parse_ref_id, print_log
from .base import ShellTask


class FetchGenomeFASTA(ShellTask):
    ref_fa_list = luigi.ListParameter()
    p = luigi.DictParameter()
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            self.ref_fa_list[0]['src'] if (
                len(self.ref_fa_list) == 1
                and (not self.ref_fa_list[0]['is_url'])
                and (str(Path(self.ref_fa_list[0]['src']).parent)
                     == self.p['ref_dir_path'])
                and self.ref_fa_list[0]['src'].endswith(('.fa', '.fasta'))
            ) else str(
                Path(self.p['ref_dir_path']).joinpath(
                    '.'.join([
                        Path(Path(Path(d['src']).name).stem).stem
                        for d in self.ref_fa_list
                    ]) + '.fa'
                )
            )
        )

    def run(self):
        fa_path = self.output().path
        run_id = parse_ref_id(ref_fa_path=fa_path)
        print_log(f'Create a reference FASTA:\t{run_id}')
        cat = self.p['cat']
        curl = self.p['curl']
        pigz = self.p['pigz']
        pbzip2 = self.p['pbzip2']
        n_cpu = self.p['n_cpu_per_worker']
        args = [
            f'{cat} --version',
            f'{curl} --version',
            f'{pigz} --version',
            f'{pbzip2} --version'
        ]
        input_files = list()
        for i, d in enumerate(self.ref_fa_list):
            s = d['src']
            if d['is_url']:
                a0 = f'{curl} -LS {s}'
            else:
                a0 = f'{cat} {s}'
                input_files.append(s)
            if s.endswith('.gz'):
                a1 = f' | {pigz} -p {n_cpu} -dc -'
            elif s.endswith('.bz2'):
                a1 = f' | {pbzip2} -p# {n_cpu} -dc -'
            else:
                a1 = ''
            r = '>' if i == 0 else '>>'
            args.append(f'set -eo pipefail && {a0}{a1} {r} {fa_path}')
        self.bash_c(
            args=args, input_files=input_files,
            output_files=fa_path, cwd=self.p['ref_dir_path'],
            run_id=run_id, log_dir_path=self.p['log_dir_path']
        )


@requires(FetchGenomeFASTA)
class CreateFASTAIndex(ShellTask):
    p = luigi.DictParameter()
    priority = 8

    def output(self):
        return luigi.LocalTarget(f'{self.input().path}.fai')

    def run(self):
        fa_path = self.input().path
        run_id = parse_ref_id(ref_fa_path=fa_path)
        print_log(f'Create a FASTA index:\t{run_id}')
        samtools = self.p['samtools']
        self.bash_c(
            args=[
                f'{samtools} 2>&1 | grep -e "Version:"',
                f'set -e && {samtools} faidx {fa_path}'
            ],
            input_files=fa_path, output_files=self.output().path,
            cwd=self.p['ref_dir_path'], run_id=run_id,
            log_dir_path=self.p['log_dir_path']
        )


@requires(FetchGenomeFASTA)
class CreateBWAIndexes(ShellTask):
    p = luigi.DictParameter()
    priority = 9

    def output(self):
        return [
            luigi.LocalTarget(f'{self.input().path}.{s}')
            for s in ['pac', 'bwt', 'ann', 'amb', 'sa']
        ]

    def run(self):
        fa_path = self.input().path
        run_id = parse_ref_id(ref_fa_path=fa_path)
        print_log(f'Create BWA indexes:\t{run_id}')
        bwa = self.p['bwa']
        self.bash_c(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'set -e && {bwa} index {fa_path}'
            ],
            input_files=fa_path,
            output_files=[o.path for o in self.output()],
            cwd=self.p['ref_dir_path'], run_id=run_id,
            log_dir_path=self.p['log_dir_path']
        )


if __name__ == '__main__':
    luigi.run()
