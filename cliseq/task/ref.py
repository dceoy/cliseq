#!/usr/bin/env python

import os
from itertools import product
from pathlib import Path

import luigi
from shoper import ShellOperator

from ..util import (open_readable_file, print_log, read_yml,
                    remove_files_if_they_exists, retrieve_url)


class FetchGenomeFASTA(luigi.Task):
    param_dict = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__ref_dict = self.param_dict['ref_fa']
        self.__ref_dir = Path(self.param_dict['ref_dir_path'])
        fa_suffixes = tuple([
            ''.join(t) for t in product(['.fa', '.fasta'], ['', '.gz', 'bz2'])
        ])
        ref_paths = self.__ref_dict.get('urls') or self.__ref_dict.get('paths')
        assert all([p.endswith(fa_suffixes) for p in ref_paths])
        self.__ref_fa_path = self.__ref_dir.joinpath(
            '.'.join([Path(Path(Path(p).name).stem).stem for p in ref_paths])
            + '.fa'
        )

    def output(self):
        return luigi.LocalTarget(self.__ref_fa_gz_path)

    def run(self):
        print_log('Create ref FASTA: {}'.format(self.__ref_fa_gz_path))
        self.__ref_dir.mkdir(exist_ok=True)
        if self.__ref_dict.get('urls'):
            urls = {
                u: self.__ref_dir.joinpath(Path(u).name)
                for u in self.__ref_dict['urls']
            }
            for u, o in urls.items():
                retrieve_url(url=u, output_path=o)
            src_paths = list(urls.values())
        else:
            urls = dict()
            src_paths = self.__ref_dict['paths']
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


class CreateFASTAIndex(luigi.Task):
    param_dict = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__samtools = self.param_dict['samtools']
        self.__ref_fa_path = Path(self.input().path).resolve()
        self.__ref_fai_path = self.__ref_fa_path + '.fai'
        self.sh = ShellOperator(
            log_txt='samtools_faidx.{}.sh.log.txt'.format(
                Path(Path(self.__ref_fa_path).name).stem
            ),
            quiet=True
        )

    def output(self):
        return luigi.LocalTarget(self.__ref_fai_path)

    def requires(self):
        return FetchGenomeFASTA()

    def run(self):
        print_log('Create a FASTA index.')
        self.sh.run(
            args=[
                f'{self.__samtools} 2>&1 | grep -e "Version:"',
                f'{self.__samtools} faidx {self.__ref_fa_path}'
            ],
            input_files=self.__ref_fa_path, output_files=self.__ref_fai_path
        )


class CreateBWAIndexes(luigi.Task):
    param_dict = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__bwa = self.param_dict['bwa']
        self.__ref_fa_path = Path(self.input().path).resolve()
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

    def requires(self):
        return FetchGenomeFASTA()

    def run(self):
        print_log('Create BWA indexes.')
        self.sh.run(
            args=[
                f'{self.__bwa} 2>&1 | grep -e "Version:"',
                f'{self.__bwa} index {self.__ref_fa_path}'
            ],
            input_files=self.__ref_fa_path, output_files=self.__ref_index_paths
        )


if __name__ == '__main__':
    param_dict = read_yml(path='task_param.yml')
    luigi.build(
        [
            CreateFASTAIndex(param_dict=param_dict),
            CreateBWAIndexes(param_dict=param_dict)
        ],
        workers=2, local_scheduler=True
    )
