#!/usr/bin/env python

from setuptools import find_packages, setup

from vcline import __version__

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='vcline',
    version=__version__,
    author='Daichi Narushima',
    author_email='dnarsil+github@gmail.com',
    description='Analytical pipeline for clinical DNA sequencing data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/dceoy/vcline.git',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['docopt', 'ftarc', 'luigi', 'psutil', 'vanqc'],
    entry_points={
        'console_scripts': ['vcline=vcline.cli.main:main']
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development'
    ],
    python_requires='>=3.6',
)
