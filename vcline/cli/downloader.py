#!/usr/bin/env python

import os
from pathlib import Path
from urllib.request import urlcleanup, urlopen, urlretrieve

import luigi
from psutil import cpu_count

from ..task.download import DownloadVCFAndExtractAF
from .util import fetch_executable, print_log, read_yml


def download_gnomad_vcf_and_extract_af(vcf_bgz_url, work_dir_path='.',
                                       max_n_cpu=None, log_level='WARNING'):
    total_n_cpu = cpu_count()
    luigi.build(
        [
            DownloadVCFAndExtractAF(
                vcf_gz_url=(
                    vcf_bgz_url or _load_default_url_dict()['gnomad_vcf']
                ),
                dest_dir_path=str(Path(work_dir_path or '.').resolve()),
                n_cpu=min(total_n_cpu, int(max_n_cpu or total_n_cpu)),
                **{c: fetch_executable(c) for c in ['curl', 'sed', 'bgzip']}
            )
        ],
        local_scheduler=True, log_level=log_level
    )


def _load_default_url_dict():
    return read_yml(
        path=str(Path(__file__).parent.joinpath('../static/urls.yml'))
    )


def download_hg(work_dir_path=None):
    print_log(f'Download genome FASTA:\thg38')
    urls = _load_default_url_dict()['ref_fa']
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
