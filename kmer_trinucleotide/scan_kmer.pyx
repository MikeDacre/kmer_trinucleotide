# -*- coding: utf-8 -*-
"""
Cython optimized kmer scanning functions
"""
import numpy as np
import pandas as pd

################################################################################
#                Wrapper Functions for Testing, Will Be Deleted                #
################################################################################

def gt(kmer, pos):
    return get_trinucleotide(kmer, pos)

def gl(kmer):
    return trinucleoties(kmer)


################################################################################
#                             Nucleotide grabbers                              #
################################################################################


cdef list trinucleoties(str kmer):
    """Return all trinucleotides from kmer as a list."""
    cdef list trios = []
    for pos in range(len(kmer)):
        our_tri = kmer[pos:pos+3]
        # Only consider complete trinucleotides
        if len(our_tri) == 3:
            trios.append((pos, our_tri))
    return trios


cdef str get_trinucleotide(str kmer, int position):
    """Get the trinucleotide from kmer at position."""
    if position > len(kmer)-3:
        return None
    return kmer[position:position+3]


################################################################################
#                            Build Kmer Dictionary                             #
################################################################################


def build_dict(list kmers):
    """Create a dictionary of kmer mean weights.

    Parameters
    ----------
    kmers : list
        List of [(kmer, weight)]

    Returns
    -------
    trinucleotides : dict
        {position : {trinucleotide: average_weight}}
    """
    cdef set tris = {
        'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC',
        'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT',
        'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC',
        'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
        'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC',
        'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT',
        'TTA', 'TTC', 'TTG', 'TTT'
    }
    cdef int kmer_range = len(kmers[0][0])-2
    cdef dict lists = {i: {k: [] for k in tris} for i in range(kmer_range)}
    cdef dict results = {i: {k: np.nan for k in tris} for i in range(kmer_range)}
    for kmer, weight in kmers:
        for pos, tri in trinucleoties(kmer):
            lists[pos][tri].append(weight)
    for pos, info in lists.items():
        for tri, weights in info.items():
            if weights:
                results[pos][tri] = np.mean(weights)
    return results


################################################################################
#                                 Get Weights                                  #
################################################################################

def get_other_weights(int pos, str tri, dict data):
    for t, w in data[pos].items():
        if t[0] == tri[0] and t[1] != tri[1] and t[2] == tri[2]:
            yield(t, w)


def get_weights(list kmers):
    """Create a dictionary of weight differences for every trinucleotide.

    Averaged per position.

    Weight of second kmer - weight of first, i.e. how does weight change when
    changing from the second trinucleotide to the first.

    Parameters
    ----------
    kmers : list
        List of [(kmer, weight)]

    Returns
    -------
    trinucleotide_weights : dict
        {trinucleotide1: {trinucleotide2: weight_change}}
        weight_change is trinucleotide2-trinucleotide1, or what the weight
        change would be if going from trinucleotide2 to trinucleotide1
    """
    cdef list nucs = list('ATGC')
    cdef dict data = build_dict(kmers)
    cdef dict trinucs = {}
    cdef dict weights = {}
    # Initialize dictionary
    for nuc1 in nucs:
        for nuc2 in nucs:
            for nuc3 in nucs:
                core = ''.join([nuc1, nuc2, nuc3])
                trinucs[core] = np.nan
                alt = {}
                for nuc_alt in [nuc for nuc in nucs if nuc != nuc2]:
                    alt[''.join([nuc1, nuc_alt, nuc3])] = np.nan
                weights[core] = alt
    # Get all mean weights
    for kmer, weight in kmers:
        for pos, tri in trinucleoties(kmer):
            for t, w in get_other_weights(pos, tri, data):
                # Primary trinucleotide weight
                orig = trinucs[tri]
                if orig is np.nan:
                    trinucs[tri] = weight
                else:
                    trinucs[tri] = np.mean([orig, weight])
                # Secondary trinucleotide weight
                orig = weights[tri][t]
                if orig is np.nan:
                    weights[tri][t] = w
                else:
                    weights[tri][t] = np.mean([orig, w])
    # Do the subtraction
    for nuc1, info in weights.items():
        for nuc2, w in info.items():
            weights[nuc1][nuc2] = w - trinucs[nuc1]
    return weights

def get_weight_matrix(list kmers):
    """Return a DataFrame of trinucleotide relative weights.

    Weight are rows minus columns.

    Averaged per position.

    Weight of second kmer - weight of first, i.e. how does weight change when
    changing from the second trinucleotide to the first.

    Parameters
    ----------
    kmers : list
        List of [(kmer, weight)]

    Returns
    -------
    DataFrame
    """
    weights = get_weights(kmers)
    df = pd.DataFrame.from_dict(weights)
