#!/usr/bin/env python

from pathlib import Path

import luigi
from luigi.util import requires

from ..cli.util import print_log
from .base import ShellTask


class FetchGenomeFASTA(ShellTask):
    log_dir_path = luigi.Parameter()
    ref_dir_path = luigi.Parameter()
    refs = luigi.ListParameter()
    ref_fa_path = luigi.Parameter()
    cat = luigi.Parameter()
    curl = luigi.Parameter()
    pigz = luigi.Parameter()
    pbzip2 = luigi.Parameter()
    n_cpu = luigi.IntParameter()
    priority = 10

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_id = parse_ref_id(ref_fa_path=self.ref_fa_path)

    def output(self):
        return luigi.LocalTarget(self.ref_fa_path)

    def run(self):
        print_log(f'Create a reference FASTA:\t{self.__ref_id}')
        self.init_bash(
            run_id=self.__ref_id, run_dir_path=self.ref_dir_path,
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
            input_files=[p for p in self.refs if not self._is_url(p)],
            output_files=self.ref_fa_path
        )

    def _geneate_commands(self):
        for i, u in enumerate(self.refs):
            if self._is_url(u):
                a0 = f'{self.curl} -LS {u}'
            else:
                a0 = f'{self.cat} ' + str(Path(u).resolve())
            if u.endswith('.gz'):
                a1 = f' | {self.pigz} -p {self.n_cpu} -dc -'
            elif u.endswith('.bz2'):
                a1 = f' | {self.pbzip2} -p# {self.n_cpu} -dc -'
            else:
                a1 = f''
            r = '>' if i == 0 else '>>'
            yield f'set -eo pipefail && {a0}{a1} {r} {self.ref_fa_path}'

    @staticmethod
    def _is_url(arg):
        return arg.startswith(('https://', 'http://'))


@requires(FetchGenomeFASTA)
class CreateFASTAIndex(ShellTask):
    log_dir_path = luigi.Parameter()
    ref_dir_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    samtools = luigi.Parameter()
    priority = 8

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_id = parse_ref_id(ref_fa_path=self.ref_fa_path)
        self.__ref_fai_path = f'{self.ref_fa_path}.fai'

    def output(self):
        return luigi.LocalTarget(self.__ref_fai_path)

    def run(self):
        print_log(f'Create a FASTA index:\t{self.__ref_id}')
        self.init_bash(
            run_id=self.__ref_id, run_dir_path=self.ref_dir_path,
            log_dir_path=self.log_dir_path
        )
        self.bash_c(
            args=[
                f'{self.samtools} 2>&1 | grep -e "Version:"',
                f'set -e && {self.samtools} faidx {self.ref_fa_path}'
            ],
            input_files=self.ref_fa_path, output_files=self.__ref_fai_path
        )


@requires(FetchGenomeFASTA)
class CreateBWAIndexes(ShellTask):
    log_dir_path = luigi.Parameter()
    ref_dir_path = luigi.Parameter()
    ref_fa_path = luigi.Parameter()
    bwa = luigi.Parameter()
    priority = 9

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_id = parse_ref_id(ref_fa_path=self.ref_fa_path)
        self.__ref_index_paths = [
            f'{self.ref_fa_path}.{s}'
            for s in ['pac', 'bwt', 'ann', 'amb', 'sa']
        ]

    def output(self):
        return luigi.LocalTarget(self.__ref_index_paths)

    def run(self):
        print_log(f'Create BWA indexes:\t{self.__ref_id}')
        self.init_bash(
            run_id=self.__ref_id, run_dir_path=self.ref_dir_path,
            log_dir_path=self.log_dir_path
        )
        self.bash_c(
            args=[
                f'{self.bwa} 2>&1 | grep -e "Version:"',
                f'set -e && {self.bwa} index {self.ref_fa_path}'
            ],
            input_files=self.ref_fa_path, output_files=self.__ref_index_paths
        )


def parse_ref_id(ref_fa_path):
    return Path(Path(ref_fa_path).name).stem


if __name__ == '__main__':
    luigi.run()
