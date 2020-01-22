#!/usr/bin/env python

import os
from pathlib import Path
from urllib.request import urlcleanup, urlopen, urlretrieve

from .util import print_log, read_yml


def download_hg(hg_ver='hg38', work_dir_path=None):
    print_log(f'Download genome FASTA:\t{hg_ver}')
    urls = read_yml(
        path=str(Path(__file__).parent.joinpath('../static/urls.yml'))
    )['ref_fa_gz'][hg_ver]
    assert all([u.endswith('.gz') for u in urls]), 'invalid gzip URLs'
    output_path = Path(work_dir_path or '.').joinpath(
        '.'.join([Path(Path(u).stem).stem for u in urls]) + '.fa.gz'
    )
    if len(urls) == 1:
        _download_file(url=urls[0], output_path=output_path)
    else:
        _download_and_merge_files(urls=urls, output_path=output_path)


def _download_file(url, output_path):
    try:
        print_log(f'Retrieve:{os.linesep}{url} => {output_path}')
        urlretrieve(url, filename=output_path)
        urlcleanup()
    except Exception as e:
        if Path(output_path).exists():
            os.remove(output_path)
        raise e


def _download_and_merge_files(urls, output_path, mode='wb'):
    try:
        with open(output_path, mode) as f:
            for u in urls:
                print_log(f'Write:{os.linesep}{u}{output_path}')
                f.write(urlopen(u).read())
    except Exception as e:
        if Path(output_path).exists():
            os.remove(output_path)
        raise e
