#!/usr/bin/env python

from setuptools import find_packages, setup

from cliseq import __version__

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='cliseq',
    version=__version__,
    author='Daichi Narushima',
    author_email='dnarsil+github@gmail.com',
    description='Analytical pipeline for clinical DNA sequencing data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='git@github.com:dceoy/cliseq.git',
    packages=find_packages(),
    include_package_data=False,
    install_requires=['docopt', 'luigi'],
    entry_points={
        'console_scripts': ['cliseq=cliseq.cli:main']
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development'
    ],
    python_requires='>=3.6',
)
