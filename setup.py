# -*- coding: utf-8 -*-
"""
Setup Script for LD_Direction
"""
import os
import codecs

from setuptools import setup

from Cython.Build import cythonize

###############################################################################
#                     Build the things we need for setup                      #
###############################################################################

# Get the long description from the README file
here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

###############################################################################
#                                Setup Options                                #
###############################################################################

setup(
    name='kmer_trinucleotide',
    version='0.1.0a1',
    description=('Lookup LD between any two lists of SNPs'),
    long_description=long_description,
    url='https://github.com/MikeDacre/kmer_trinucleotide',
    author='Michael Dacre',
    author_email='mike.dacre@gmail.com',
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Operating System :: Linux',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    keywords='genetics',

    install_requires=['cython'],
    packages=['kmer_trinucleotide'],
    #  entry_points={
        #  'console_scripts': [
            #  'ldpair = LD_Direction.LDpair:main',
            #  'ldlists = LD_Direction.compare_snp_lists:main',
        #  ],
    #  },
    ext_modules = cythonize('kmer_trinucleotide/scan_kmer.pyx'),
)
