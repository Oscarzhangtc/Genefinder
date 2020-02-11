# -*- coding: utf-8 -*-
"""
GeneFinder Week 2

@author: Oscar Zhang

"""

import random
import doctest
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

"""
# uncomment to run without loading the data in terminal
 dna = load_seq("./data/X73525.fa")
"""

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

    additional unit tests:  test out all possibilities
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
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
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_list_complement=[]
        # set list to a blank array
    for letter in dna:
        complement = get_complement(letter)
        # call the get_complement_of_ORf element functoin prior
        reverse_list_complement.insert(0, complement)
        reverse_complement = ''.join(reverse_list_complement)

    return reverse_complement
    # TODO: implement this
    pass

# doctest.run_docstring_examples(get_reverse_complement, globals(), verbose=False)

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

    additional unit test: this additional unit test ensure the while loop stops at the right stop codon when there are multiple stop codons
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGTGAATGA")
    'ATG'
    """
    i = 0
    stop = len(dna)
    while i <= (stop/3):
        # here we basically divide the entire dna into sets of 3s
        # we can run the loop number of set times
        # e.g if there are 4 sets of characters(total length of 12 characters), we would run the loop 4 times
        if (dna[ 3 * i : 3 * i + 3]== 'TAG' or dna[ 3 * i : 3 * i + 3]=='TGA' or dna[3*i:3*i+3]=='TAA'):
            stop = i*3
            # run through evey set of 3, and record the index when a set = stop codon
            break

        else:
            i+=1
    return dna[:stop]
    # return the dna string up to the stopping point

    # TODO: implement this
    pass

# doctest.run_docstring_examples(rest_of_ORF, globals(), verbose=False)

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    additional unit test: ensures the returned frame always starts at the 1st start codon
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGATGCATGAATGTAGATAGATGTGCCC")
    ['ATGATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # initialize the start points, and frame
    one_frame = []
    start_codon = 'ATG'
    start_point = -1
    while len(dna)/3 >= 1:
        # ensure length of dna is larger than 3
        i = 0
        while i <= len(dna)/3:
            if (dna[i * 3: i * 3 +3] == start_codon):
                start_point = 3*i
                break
            else:
                i += 1
            # not breaking until we find the start codon
        if (start_point!=-1):
            # continues if we find start codon
            dna = dna[start_point:]
            # returns the remainder of the dna from its start codon
        else:
            return one_frame
        orf = rest_of_ORF(dna)
        # orf = dna starting with a start_codon up to its stop codon(not included)
        cut = len(dna)-len(orf)-3
        # how much is cut off
        start_point = len(dna)- cut
        # reset start point to equal to the place its cut off
        dna = dna[start_point:]
        # reset dna
        one_frame.append(orf)
        # append orf to one_frame list
        # the process repeats as long the remaining dna is > 3 characters
    one_frame = list(filter(None, one_frame))
    return one_frame
    # TODO: implement this
    pass


# doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=False)

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    additional unit test: ensures function is able to consistently find the start codons
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("TAGATGCATGAATGTAGTAG")
    ['ATGCATGAATGTAGTAG', 'ATGAATGTAGTAG', 'ATG']
    """
    orf_list = []
    for a in range(0 , 3):
        # runs 3 times to find all 3 possible frams
        dna1= dna[a: ]
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
    no additional unit test: the function relies only functions previously tested and verified
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    orf_list = find_all_ORFs(dna)
    reverse_complement = get_reverse_complement(dna)
    reverse_orf_list = find_all_ORFs(reverse_complement)
    orf_list.extend(reverse_orf_list)
    # combine the resulting orf_lists in to a single list/array

    return orf_list
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """

    orf_list = find_all_ORFs_both_strands(dna)
    longest_orf_list = max(orf_list, key = len)
    # performs max() on all elements in orf_list
    # max() finds the longest element
    return longest_orf_list
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    output = []
    for i in range (num_trials):
        test_dna = shuffle_string(dna)
        # random shuffles dna
        longest_ORF(dna)
        output.extend(dna)
        # combines the resuliting long_ORF(s) in to a single list
    longest_ORF_noncoding = max(output,key = len)
    # apply max() to find longest_ORF_noncoding
    return len(longest_ORF_noncoding)
    # TODO: implement this
    pass

# doctest.run_docstring_examples(longest_ORF_noncoding, globals(), verbose=False)
# error-prone code, added doctest to ensure correctness
def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        no additional unit test is required the function mainly relies on the aa_table method and all logic are previously verified
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    amino_acid_list = []
    stop = int(len(dna)/3)
    for i in range (stop):
        piece = dna[i*3:i*3+3]
        amino_acid = aa_table[piece]
        # use aa_table funtion imported from "amino_acid" folder and apply to the set of 3
        amino_acid_list.extend(amino_acid)
    amino_acid_list = ''.join(amino_acid_list)
    # concatenates into string
    return amino_acid_list
    # TODO: implement this
    pass

# doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=False)

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    sequence_list=[]
    threshold = longest_ORF_noncoding(dna, 1500)
    # shuffling 1500 times and find the longest orf
    all_orf = find_all_ORFs_both_strands(dna)
    # the
    for orf in all_orf:
        if (len(orf) >= threshold):
            #length of the return orfs is larger than the threshhold or the length of the longest ORF_noncoding
            sequence_list.extend(orf)
    sequence_list = ''.join(sequence_list)
    print(coding_strand_to_AA(sequence_list))
    return coding_strand_to_AA(sequence_list)
    # computes the protein

    # TODO: implement this
    pass

"""
# uncomment run to without loading the data in terminal
gene_finder(dna)
"""

if __name__ == "__main__":
    import doctest
    doctest.testmod()
