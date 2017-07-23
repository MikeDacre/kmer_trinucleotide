# -*- coding: utf-8 -*-
"""
Compare trinucleotide effects for any list of kmers with wieghts.

Info
----
Author: Michael D Dacre, mike.dacre@gmail.com
Organization: Stanford University
License: MIT License, property of Stanford, use as you wish
Version: 0.1.0a1
Created: 2017-09-19 16:07
Last modified: 2017-07-19 16:53

Summary
-------

This algorithm will loop through any file of:

+------+--------+
| kmer | weight |
+------+--------+

For each kmer, it will loop through all available trinucleotides, and then
compare those to **every other kmer at that position**. For each trinucleotide
it will ask how, on *average* does the **central nucleotide** affect the score.
In other words, changing the central nucleotide in any trinucleotide changes
the weight by how much?

This is currently under development, as in its current form it will take several
years to run.
"""
__version__ = '0.1.0a1'

from . import scan_kmer
