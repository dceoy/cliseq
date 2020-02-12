#!/usr/bin/env python

import os
from pathlib import Path
from urllib.request import urlcleanup, urlopen, urlretrieve

from .util import load_default_url_dict, print_log


def download_hg(work_dir_path=None):
    print_log(f'Download genome FASTA:\thg38')
    urls = load_default_url_dict()['ref_fa']
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
