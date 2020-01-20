#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .bash import BashTask


class FetchGenomeFASTA(BashTask):
    log_dir_path = luigi.Parameter()
    ref_dir_path = luigi.Parameter()
    ref_fa = luigi.DictParameter()
    ref_id = luigi.Parameter()
    cat = luigi.Parameter()
    curl = luigi.Parameter()
    pigz = luigi.Parameter()
    pbzip2 = luigi.Parameter()
    n_cpu = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_fa_path = generate_ref_fa_path(
            ref_id=self.ref_id, ref_dir_path=self.ref_dir_path
        )

    def output(self):
        return luigi.LocalTarget(self.__ref_fa_path)

    def run(self):
        print_log(f'Create a reference FASTA:\t{self.ref_id}')
        self.init_bash(
            log_name=self.ref_id, run_dir_path=self.ref_dir_path,
            log_dir_path=self.log_dir_path
        )
        self.bash_c(
            args=[
                f'{self.cat} --version',
                f'{self.curl} --version',
                f'{self.pigz} --version',
                f'{self.pbzip2} --version',
                *self._geneate_commands()
            ],
            input_files=[p for p in self.ref_fa if not self._is_url(p)],
            output_files=self.__ref_fa_path
        )

    def _geneate_commands(self):
        for i, u in enumerate(self.ref_fa):
            cat = f'{self.curl} -LS' if self._is_url(u) else f'{self.cat}'
            if u.endswith('.gz'):
                a = f'{cat} {u} | {self.pigz} -p {self.n_cpu} -dc -'
            elif u.endswith('.bz2'):
                a = f'{cat} {u} | {self.pbzip2} -p# {self.n_cpu} -dc -'
            else:
                a = f'{cat} {u}'
            r = '>' if i == 0 else '>>'
            yield f'set -eo pipefail && {a} {r} {self.__ref_fa_path}'

    @staticmethod
    def _is_url(arg):
        return arg.startswith(('https://', 'http://'))


@requires(FetchGenomeFASTA)
class CreateFASTAIndex(BashTask):
    log_dir_path = luigi.Parameter()
    ref_dir_path = luigi.Parameter()
    ref_id = luigi.Parameter()
    samtools = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_fai_path = f'{self.input().path}.fai'

    def output(self):
        return luigi.LocalTarget(self.__ref_fai_path)

    def run(self):
        print_log(f'Create a FASTA index:\t{self.ref_id}')
        self.init_bash(
            log_name=self.ref_id, run_dir_path=self.ref_dir_path,
            log_dir_path=self.log_dir_path
        )
        self.bash_c(
            args=[
                f'{self.samtools} 2>&1 | grep -e "Version:"',
                f'set -e && {self.samtools} faidx {self.input().path}'
            ],
            input_files=self.input().path, output_files=self.__ref_fai_path
        )


@requires(FetchGenomeFASTA)
class CreateBWAIndexes(BashTask):
    log_dir_path = luigi.Parameter()
    ref_dir_path = luigi.Parameter()
    ref_id = luigi.Parameter()
    bwa = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_index_paths = [
            f'{self.input().path}{s}'
            for s in ['.pac', '.bwt', '.ann', '.amb', '.sa']
        ]

    def output(self):
        return luigi.LocalTarget(self.__ref_index_paths)

    def run(self):
        print_log(f'Create BWA indexes:\t{self.ref_id}')
        self.init_bash(
            log_name=self.ref_id, run_dir_path=self.ref_dir_path,
            log_dir_path=self.log_dir_path
        )
        self.bash_c(
            args=[
                f'{self.bwa} 2>&1 | grep -e "Version:"',
                f'set -e && {self.bwa} index {self.input().path}'
            ],
            input_files=self.input().path, output_files=self.__ref_index_paths
        )


def generate_ref_fa_path(ref_id, ref_dir_path):
    return Path(ref_dir_path).joinpath(f'{ref_id}.fa')


if __name__ == '__main__':
    luigi.run()
