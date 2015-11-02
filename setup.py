#!/usr/bin/env python

from setuptools import setup

setup(
    name='thepipeline',
    version='0.1',
    author='Khalid Mahmood',
    author_email='khalid.mahmood@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['thepipepine = src.main:main']
    },
    url='https://github.com/khalidm/thepipepine',
    license='LICENSE.txt',
    description='thepipepine is a bioinformatics pipeline to filter and \
     annotate a vcf file.',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
