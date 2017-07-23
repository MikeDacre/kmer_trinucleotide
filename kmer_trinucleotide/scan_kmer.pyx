# -*- coding: utf-8 -*-
"""
Cython optimized kmer scanning functions
"""

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
