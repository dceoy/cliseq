#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .base import ShellTask


class FetchGenomeFASTA(ShellTask):
    params = luigi.DictParameter()
    priority = 10

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ref_src_paths = [
            (p if _is_url(p) else str(Path(p).resolve()))
            for p in self.params['refs']
        ]
        self.ref_fa_path = self.params['ref_fa_path']

    def output(self):
        return luigi.LocalTarget(self.ref_fa_path)

    def run(self):
        ref_fa_path = self.params['ref_fa_path']
        ref_id = parse_ref_id(ref_fa_path=ref_fa_path)
        cat = self.params['cat']
        curl = self.params['curl']
        pigz = self.params['pigz']
        pbzip2 = self.params['pbzip2']
        n_cpu = self.params['n_cpu']
        print_log(f'Create a reference FASTA:\t{ref_id}')
        self.init_bash(
            run_id=ref_id, run_dir_path=self.params['ref_dir_path'],
            log_dir_path=self.params['log_dir_path']
        )
        cmds = list()
        for i, u in enumerate(self.ref_src_paths):
            a0 = f'{curl} -LS {u}' if _is_url(u) else f'{cat} {u}'
            if u.endswith('.gz'):
                a1 = f' | {pigz} -p {n_cpu} -dc -'
            elif u.endswith('.bz2'):
                a1 = f' | {pbzip2} -p# {n_cpu} -dc -'
            else:
                a1 = f''
            r = '>' if i == 0 else '>>'
            cmds.append(f'set -eo pipefail && {a0}{a1} {r} {ref_fa_path}')
        self.bash_c(
            args=[
                f'{cat} --version', f'{curl} --version', f'{pigz} --version',
                f'{pbzip2} --version', *cmds
            ],
            input_files=[p for p in self.ref_src_paths if not _is_url(p)],
            output_files=self.ref_fa_path
        )


@requires(FetchGenomeFASTA)
class CreateFASTAIndex(ShellTask):
    params = luigi.DictParameter()
    priority = 8

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ref_fa_path = self.params['ref_fa_path']
        self.ref_fai_path = f'{self.ref_fa_path}.fai'

    def output(self):
        return luigi.LocalTarget(self.ref_fai_path)

    def run(self):
        ref_id = parse_ref_id(ref_fa_path=self.ref_fa_path)
        samtools = self.params['samtools']
        print_log(f'Create a FASTA index:\t{ref_id}')
        self.init_bash(
            run_id=ref_id, run_dir_path=self.params['ref_dir_path'],
            log_dir_path=self.params['log_dir_path']
        )
        self.bash_c(
            args=[
                f'{samtools} 2>&1 | grep -e "Version:"',
                f'set -e && {samtools} faidx {self.ref_fa_path}'
            ],
            input_files=self.ref_fa_path, output_files=self.ref_fai_path
        )


@requires(FetchGenomeFASTA)
class CreateBWAIndexes(ShellTask):
    params = luigi.DictParameter()
    priority = 9

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ref_fa_path = self.params['ref_fa_path']
        self.ref_index_paths = [
            f'{self.ref_fa_path}.{s}'
            for s in ['pac', 'bwt', 'ann', 'amb', 'sa']
        ]

    def output(self):
        return luigi.LocalTarget(self.ref_index_paths)

    def run(self):
        ref_id = parse_ref_id(ref_fa_path=self.ref_fa_path)
        bwa = self.params['bwa']
        print_log(f'Create BWA indexes:\t{ref_id}')
        self.init_bash(
            run_id=ref_id, run_dir_path=self.params['ref_dir_path'],
            log_dir_path=self.params['log_dir_path']
        )
        self.bash_c(
            args=[
                f'{bwa} 2>&1 | grep -e "Version:"',
                f'set -e && {bwa} index {self.ref_fa_path}'
            ],
            input_files=self.ref_fa_path, output_files=self.ref_index_paths
        )


def _is_url(arg):
    return arg.startswith(('https://', 'http://'))


def parse_ref_id(ref_fa_path):
    return Path(Path(ref_fa_path).name).stem


if __name__ == '__main__':
    luigi.run()
