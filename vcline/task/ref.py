#!/usr/bin/env python

import os
from itertools import product
from pathlib import Path

import luigi
from luigi.util import requires
from shoper.shelloperator import ShellOperator

from ..cli.util import (download_file, open_readable_file, print_log,
                        remove_files_if_they_exists)


class FetchGenomeFASTA(luigi.Task):
    ref_fa = luigi.DictParameter()
    ref_dir_path = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_dir = Path(self.ref_dir_path)
        fa_suffixes = tuple([
            ''.join(t) for t in product(['.fa', '.fasta'], ['', '.gz', 'bz2'])
        ])
        ref_paths = self.ref_fa.get('urls') or self.ref_fa.get('paths')
        assert all([p.endswith(fa_suffixes) for p in ref_paths])
        self.__ref_fa_path = self.__ref_dir.joinpath(
            '.'.join([Path(Path(Path(p).name).stem).stem for p in ref_paths])
            + '.fa'
        )

    def output(self):
        return luigi.LocalTarget(self.__ref_fa_path)

    def run(self):
        print_log('Create ref FASTA: {}'.format(self.__ref_fa_path))
        self.__ref_dir.mkdir(exist_ok=True)
        if self.ref_fa.get('urls'):
            urls = {
                u: str(self.__ref_dir.joinpath(Path(u).name))
                for u in self.ref_fa['urls']
            }
            for u, o in urls.items():
                download_file(url=u, output_path=o)
            src_paths = list(urls.values())
        else:
            urls = dict()
            src_paths = self.ref_fa['paths']
        try:
            with open(self.__ref_fa_path, 'w') as fw:
                for p in src_paths:
                    print_log(
                        'Write:{0}{1} => {2}'.format(
                            os.linesep, p, self.__ref_fa_path
                        )
                    )
                    with open_readable_file(path=p) as fr:
                        for line in fr:
                            fw.write(line)
        except Exception as e:
            if Path(self.__ref_fa_path).exists():
                os.remove(self.__ref_fa_path)
            raise e
        else:
            print_log('Save: {}'.format(self.__ref_fa_path))
        finally:
            remove_files_if_they_exists(*urls.values())


@requires(FetchGenomeFASTA)
class CreateFASTAIndex(luigi.Task):
    samtools = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_fa_path = str(Path(self.input().path).resolve())
        self.__ref_fai_path = self.__ref_fa_path + '.fai'
        self.sh = ShellOperator(
            log_txt='samtools_faidx.{}.sh.log.txt'.format(
                Path(Path(self.__ref_fa_path).name).stem
            ),
            quiet=True
        )

    def output(self):
        return luigi.LocalTarget(self.__ref_fai_path)

    def run(self):
        print_log('Create a FASTA index.')
        self.sh.run(
            args=[
                f'{self.samtools} 2>&1 | grep -e "Version:"',
                f'{self.samtools} faidx {self.__ref_fa_path}'
            ],
            input_files=self.__ref_fa_path, output_files=self.__ref_fai_path
        )


@requires(FetchGenomeFASTA)
class CreateBWAIndexes(luigi.Task):
    bwa = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_fa_path = str(Path(self.input().path).resolve())
        self.__ref_index_paths = [
            (self.__ref_fa_path + s)
            for s in ['.pac', '.bwt', '.ann', '.amb', '.sa']
        ]
        self.sh = ShellOperator(
            log_txt='bwa_index.{}.sh.log.txt'.format(
                Path(Path(self.__ref_fa_path).name).stem
            ),
            quiet=True
        )

    def output(self):
        return luigi.LocalTarget(self.__ref_index_paths)

    def run(self):
        print_log('Create BWA indexes.')
        self.sh.run(
            args=[
                f'{self.bwa} 2>&1 | grep -e "Version:"',
                f'{self.bwa} index {self.__ref_fa_path}'
            ],
            input_files=self.__ref_fa_path, output_files=self.__ref_index_paths
        )


@requires(CreateFASTAIndex)
@requires(CreateBWAIndexes)
class PrepareReferences(luigi.Task):
    pass


if __name__ == '__main__':
    luigi.run()
