#!/usr/bin/env python

import os
from itertools import product
from pathlib import Path

import luigi
from luigi.util import requires
from shoper.shelloperator import ShellOperator

from ..cli.util import print_log


class FetchGenomeFASTA(luigi.Task):
    ref_fa = luigi.DictParameter()
    ref_dir_path = luigi.Parameter()
    pigz = luigi.Parameter()
    pbzip2 = luigi.Parameter()
    n_cpu = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__log_dir = Path(self.log_dir_path)
        self.__ref_dir = Path(self.ref_dir_path)
        fa_suffixes = tuple([
            ''.join(t) for t in product(['.fa', '.fasta'], ['', '.gz', 'bz2'])
        ])
        assert all([p.endswith(fa_suffixes) for p in self.ref_fa])
        self.__ref_fa_path = self.__ref_dir.joinpath(
            '.'.join([Path(Path(Path(p).name).stem).stem for p in self.ref_fa])
            + '.fa'
        )

    def output(self):
        return luigi.LocalTarget(self.__ref_fa_path)

    def run(self):
        self.__log_dir.mkdir(exist_ok=True)
        self.__ref_dir.mkdir(exist_ok=True)
        print_log(os.linesep.join(['Create a reference FASTA:'] + self.ref_fa))
        sh = ShellOperator(
            log_txt=str(
                self.__log_dir.joinpath(
                    'cat_ref.{}.sh.log.txt'.format(
                        Path(Path(self.__ref_fa_path).name).stem
                    )
                )
            ),
            quiet=True
        )
        pigz = f'self.pigz -p {self.n_cpu}'
        pbzip2 = f'self.pbzip2 -p# {self.n_cpu}'
        ref_urls = {
            u: str(self.__ref_dir.joinpath(Path(u).name))
            for u in self.ref_fa if u.startswith(('https://', 'http://'))
        }
        args = list()
        for i, u in enumerate(self.ref_fa):
            cat = 'curl -LS' if u in ref_urls else 'cat'
            r = '>' if i == 0 else '>>'
            if u.endswith('.gz'):
                a = f'{cat} {u} | {pigz} -dc - {r} {self.__ref_fa_path}'
            elif u.endswith('.bz2'):
                a = f'{cat} {u} | {pbzip2} -dc - {r} {self.__ref_fa_path}'
            else:
                a = f'{cat} {u} {r} {self.__ref_fa_path}'
            args.append(a)
        sh.run(
            args=args,
            input_files=[p for p in self.ref_fa if p not in ref_urls],
            output_files=self.__ref_fa_path,
            cwd=str(Path(self.__ref_fa_path).parent)
        )


@requires(FetchGenomeFASTA)
class CreateFASTAIndex(luigi.Task):
    log_dir_path = luigi.Parameter()
    samtools = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_fa_path = str(Path(self.input().path).resolve())
        self.__ref_fai_path = self.__ref_fa_path + '.fai'

    def output(self):
        return luigi.LocalTarget(self.__ref_fai_path)

    def run(self):
        print_log('Create a FASTA index:\t{}'.format(self.__ref_fa_path))
        sh = ShellOperator(
            log_txt=str(
                Path(self.log_dir_path).joinpath(
                    'samtools_faidx.{}.sh.log.txt'.format(
                        Path(Path(self.__ref_fa_path).name).stem
                    )
                )
            ),
            quiet=True
        )
        sh.run(
            args=[
                f'{self.samtools} 2>&1 | grep -e "Version:"',
                f'{self.samtools} faidx {self.__ref_fa_path}'
            ],
            input_files=self.__ref_fa_path,
            output_files=self.__ref_fai_path,
            cwd=str(Path(self.__ref_fa_path).parent)
        )


@requires(FetchGenomeFASTA)
class CreateBWAIndexes(luigi.Task):
    log_dir_path = luigi.Parameter()
    bwa = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_fa_path = str(Path(self.input().path).resolve())
        self.__ref_index_paths = [
            (self.__ref_fa_path + s)
            for s in ['.pac', '.bwt', '.ann', '.amb', '.sa']
        ]

    def output(self):
        return luigi.LocalTarget(self.__ref_index_paths)

    def run(self):
        print_log('Create BWA indexes:\t{}'.format(self.__ref_fa_path))
        sh = ShellOperator(
            log_txt=str(
                Path(self.log_dir_path).joinpath(
                    'bwa_index.{}.sh.log.txt'.format(
                        Path(Path(self.__ref_fa_path).name).stem
                    )
                )
            ),
            quiet=True
        )
        sh.run(
            args=[
                f'{self.bwa} 2>&1 | grep -e "Version:"',
                f'{self.bwa} index {self.__ref_fa_path}'
            ],
            input_files=self.__ref_fa_path,
            output_files=self.__ref_index_paths,
            cwd=str(Path(self.__ref_fa_path).parent)
        )


if __name__ == '__main__':
    luigi.run()
