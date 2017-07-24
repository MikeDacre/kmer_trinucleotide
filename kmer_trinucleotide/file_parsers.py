# -*- coding: utf-8 -*-
"""Parse different files directly into a trinucleotide matrix.

Functions
---------
parse_kmer_file(filename, kmer_col, weight_col, sep='\\t')
    Create a dictionary of trinucleotide relative weights from a columnar
    file.
"""
import bz2 as _bz2
import gzip as _gzip

from numpy import float128 as flt

from . import scan_kmer as _scan


def parse_kmer_file(filename, kmer_col, weight_col, header=True, sep='\t'):
    """Parse a columnar file to build trinucleotide matrix.
import bz2 as _bz2
import gzip as _gzip

    Parameters
    ----------
    filename : str
        Path to the file to parse
    kmer_col : int
        0-based integer column number for the kmer column
    weight_col : int
        0-based integer column number for the weight column
    header : bool
        If True, ignore first line
    sep : str
        Column to split the file on

    Returns
    -------
    trinucleotide_weights : dict
        {trinucleotide1: {trinucleotide2: weight_change}}
        weight_change is trinucleotide2-trinucleotide1, or what the weight
        change would be if going from trinucleotide2 to trinucleotide1
    """
    kmers = []
    with _open_zipped(filename) as fin:
        if header:
            fin.readline()
        for l in fin:
            f = l.strip().split(sep)
            kmers.append((f[kmer_col], flt(f[weight_col])))
    return _scan.get_weights(kmers)


def _open_zipped(infile, mode='r'):
    """Return file handle of file regardless of compressed or not.

    Also returns already opened files unchanged, text mode automatic for
    compatibility with python2.
    """
    # Return already open files
    if hasattr(infile, 'write'):
        return infile
    # Make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'
    if not isinstance(infile, str):
        raise ValueError("I cannot open a filename that isn't a string.")
    if infile.endswith('.gz'):
        return _gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(_bz2, 'open'):
            return _bz2.open(infile, mode)
        else:
            return _bz2.BZ2File(infile, mode)
    return open(infile, mode)

