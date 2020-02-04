# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Oscar Zhang

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if (nucleotide == 'A'):
        complement = 'T'
    if (nucleotide == 'T'):
        complement = 'A'
    if (nucleotide == 'C'):
        complement = 'G'
    if (nucleotide == 'G'):
        complement = 'C'
    return complement
    # TODO: implement this
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverserest_of_ORf complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'rest_of_ORf
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_list_complement=[]
        # set list to a blank array
    for letter in dna:

        complement = get_complement(letter)
        # call the get_comperest_of_ORflement functoin prior
        reverse_list_complement.insert(0,complement)
        reverse_complement = ''.join(list_reverse_complement)
    return reverse_complement
    # TODO: implement this
    pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    i=0
    stop = len(dna)
    while i <= (stop/3):
        if (dna[ 3 * i : 3 * i + 3]== 'TAG' or dna[ 3 * i : 3 * i + 3]=='TGA' or dna[3*i:3*i+3]=='TAA'):
            stop = i*3
            break

        else:
            i+=1
    return dna[:stop]

    # TODO: implement this
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    one_frame = []
    start_codon = 'ATG'
    start_point = -1
    while len(dna) / 3 >= 1:
        i = 0
        while i < = len (dna) / 3:
            if (dna[i * 3 : i * 3 + 3] == start_codon):
                start_point = 3 * i
                break
            else:
                i += 1
        if (start_point != -1):
            dna = dna[start_point:]
        else:
            return one_frame
        orf = rest_of_ORF(dna)
        rest_len = len(dna) - len(orf) -3
        start_point = len(dna) - rest_len
        dna = dna[start_point:]
        one_frame.append(orf)
    one_frame = list(filter(None, one_frame))
    return one_frame
    # TODO: implement this
    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    orf_list = []
    for a in range(0,3):
        dna1=dna[a:]
        orf1 = find_all_ORFs_oneframe(dna1)
        orf_list.extend(orf1)

    return orf_list
    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # initialize the start points, and frame
    one_frame = []
    start_codon = 'ATG'
    start_point = -1
    while len(dna) /3 > 1:
        i = 0
        while i <= len(dna)/3:
            if (dna[i * 3: i * 3 +3] == start_codon):
                start_point = 3*i
                break
            else:
                i += 1
        if (start_point!=-1):
            dna = dna[start_point:]
        else:
            return one_frame
        orf = rest_of_ORF(dna)
        rest_len = len(dna)-len(orf)-3
        start_point = len(dna)-rest_len
        dna = dna[start_point:]
        one_frame.append(orf)
    one_frame = list(filter(None, one_frame))
    return one_frame
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
