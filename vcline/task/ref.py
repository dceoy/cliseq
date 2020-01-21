#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import is_url, parse_ref_id, print_log
from .base import ShellTask


class FetchGenomeFASTA(ShellTask):
    p = luigi.DictParameter()
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            str(
                Path(self.p['ref_dir_path']).joinpath(
                    '.'.join([
                        Path(Path(Path(p).name).stem).stem
                        for p in self.p['ref_fa']
                    ]) + '.fa'
                )
            )
        )

    def run(self):
        ref_fa_path = self.output().path
        ref_id = parse_ref_id(ref_fa_path=ref_fa_path)
        print_log(f'Create a reference FASTA:\t{ref_id}')
        cat = self.p['cat']
        curl = self.p['curl']
        pigz = self.p['pigz']
        pbzip2 = self.p['pbzip2']
        n_cpu = self.p['n_cpu']
        args = [
            f'{cat} --version', f'{curl} --version', f'{pigz} --version',
            f'{pbzip2} --version'
        ]
        for i, u in enumerate(self.p['ref_fa']):
            if is_url(u):
                a0 = f'{curl} -LS {u}'
            else:
                a0 = f'{cat} ' + str(Path(u).resolve())
            if u.endswith('.gz'):
                a1 = f' | {pigz} -p {n_cpu} -dc -'
            elif u.endswith('.bz2'):
                a1 = f' | {pbzip2} -p# {n_cpu} -dc -'
            else:
                a1 = f''
            r = '> ' if i == 0 else '>> '
            args.append(f'set -eo pipefail && {a0}{a1} {r} {ref_fa_path}')
        self.bash_c(
            args=args,
            input_files=[p for p in self.p['ref_fa'] if not is_url(p)],
            output_files=ref_fa_path, cwd=self.p['ref_dir_path'],
            run_id=ref_id, log_dir_path=self.p['log_dir_path']
        )


@requires(FetchGenomeFASTA)
class CreateFASTAIndex(ShellTask):
    p = luigi.DictParameter()
    priority = 8

    def output(self):
        return luigi.LocalTarget(f'{self.input().path}.fai')

    def run(self):
        ref_fa_path = self.input().path
        ref_id = parse_ref_id(ref_fa_path=ref_fa_path)
        print_log(f'Create a FASTA index:\t{ref_id}')
        samtools = self.p['samtools']
        self.bash_c(
            args=[
                f'{samtools} 2>&1 | grep -e "Version:"',
                f'set -e && {samtools} faidx {ref_fa_path}'
            ],
            input_files=ref_fa_path, output_files=self.output().path,
            cwd=self.p['ref_dir_path'], run_id=ref_id,
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
        ref_fa_path = self.input().path
        ref_id = parse_ref_id(ref_fa_path=ref_fa_path)
        print_log(f'Create BWA indexes:\t{ref_id}')
        bwa = self.p['bwa']
        self.bash_c(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'set -e && {bwa} index {ref_fa_path}'
            ],
            input_files=ref_fa_path,
            output_files=[o.path for o in self.output()],
            cwd=self.p['ref_dir_path'], run_id=ref_id,
            log_dir_path=self.p['log_dir_path']
        )


if __name__ == '__main__':
    luigi.run()
