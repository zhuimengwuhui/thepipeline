#!/usr/bin/env python

from setuptools import setup

setup(
    name='dog_pipeline',
    version='0.1',
    author='Khalid Mahmood',
    author_email='khalid.mahmood@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['dog_pipeline = src.main:main']
    },
    url='https://github.com/khalidm/thepipepine',
    license='LICENSE.txt',
    description='dog_pipeline is a bioinformatics pipeline to call variants in dog exome data.',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
