#!/usr/bin/env python

import logging
import os
import re
import subprocess
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader


def print_log(message):
    logger = logging.getLogger(__name__)
    logger.debug(message)
    print(f'>>\t{message}', flush=True)


def parse_fq_id(fq_path):
    return (
        re.sub(
            r'([\._]R[12][\._]|\.fq\.|\.fastq\.).+$', '',
            Path(fq_path).name
        ) or Path(Path(fq_path).stem).stem
    )


def create_matched_id(foreground_name, background_name):
    fragments = [
        Path(n).stem.split('.') for n in [foreground_name, background_name]
    ]
    if fragments[0][-1] != fragments[1][-1]:
        somatic_id = '.'.join(['.'.join(f) for f in fragments])
    else:
        n_common = 0
        for i in range(1, min([len(f) for f in fragments])):
            if fragments[0][-i] == fragments[1][-i]:
                n_common += 1
            else:
                break
        somatic_id = '.'.join(fragments[0][:-n_common] + fragments[1])
    return somatic_id


def fetch_executable(cmd):
    executables = [
        cp for cp in [
            str(Path(p).joinpath(cmd))
            for p in os.environ['PATH'].split(os.pathsep)
        ] if os.access(cp, os.X_OK)
    ]
    if executables:
        return executables[0]
    else:
        raise RuntimeError(f'command not found: {cmd}')


def read_yml(path):
    with open(path, 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    return d


def render_template(template, data, output_path):
    print_log(
        ('Overwrite' if Path(output_path).exists() else 'Render')
        + f' a file:\t{output_path}'
    )
    with open(output_path, 'w') as f:
        f.write(
            Environment(
                loader=FileSystemLoader(
                    str(Path(__file__).parent.joinpath('../template')),
                    encoding='utf8'
                )
            ).get_template(template).render(data) + os.linesep
        )


def _run_and_parse_subprocess(args, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, **kwargs):
    logger = logging.getLogger(__name__)
    logger.debug(f'args: {args}')
    with subprocess.Popen(args=args, stdout=stdout, stderr=stderr,
                          **kwargs) as p:
        for line in p.stdout:
            yield line.decode('utf-8')
        outs, errs = p.communicate()
        if p.returncode != 0:
            logger.error(
                'STDERR from subprocess `{0}`:{1}{2}'.format(
                    p.args, os.linesep, errs.decode('utf-8')
                )
            )
            raise subprocess.CalledProcessError(
                returncode=p.returncode, cmd=p.args, output=outs,
                stderr=errs
            )


def curl_and_write_af_only_vcf_bgz(vcf_bgz_url, vcf_bgz_path, curl, bgzip,
                                   n_cpu=1, shell=True, executable='/bin/bash',
                                   **kwargs):
    logger = logging.getLogger(__name__)
    search_regex = re.compile('\tPASS\t.*[\t;]AF=[^;]')
    sub_regexes = [
        re.compile('(\t[^\t]*;|\t)(AF=[0-9]*\\.[e0-9+-]*)[^\t\r\n]*'), '\t\\2'
    ]
    args0 = f'{curl} -LS {vcf_bgz_url} | {bgzip} -@ {n_cpu} -dc -'
    args1 = f'{bgzip} -@ {n_cpu} -c > {vcf_bgz_path}'
    logger.info(f'`{args0}` -> (process) -> `{args1}`')
    popen_kwargs = {'shell': shell, 'executable': executable, **kwargs}
    logger.debug(f'popen_kwargs:\t{popen_kwargs}')
    p1 = subprocess.Popen(
        args=args1, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, **popen_kwargs
    )
    for s in _run_and_parse_subprocess(args=args0, **popen_kwargs):
        if search_regex.search(s):
            p1.stdin.write(re.sub(*sub_regexes, s).encode('utf-8'))
        elif s.startswith('#'):
            p1.stdin.write(s.encode('utf-8'))
    p1.stdin.flush()
    p1.stdin.close()
    p1.wait()
    if p1.returncode != 0:
        logger.error(
            f'STDERR from subprocess `{p1.args}`:'
            + os.linesep + p1.stderr.decode('utf-8')
        )
        raise subprocess.CalledProcessError(
            returncode=p1.returncode, cmd=p1.args, output=p1.stdout,
            stderr=p1.stderr
        )
    elif not Path(vcf_bgz_path).is_file():
        raise FileNotFoundError(vcf_bgz_path)
