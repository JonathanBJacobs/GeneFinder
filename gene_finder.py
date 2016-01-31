# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Jonathan Jacobs

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

        The two added unit tests test the other two possible inputs that the
        function will handle. That means every viable input is tested

        The last unit test makes sure that the funciton doesnt have a fit if it
        gets a character that it shouldn't
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    >>> get_complement('B')

    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'
    else:
        return None


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        I don't believe any other unit tests are necessary because the two
        following unit tests test all four nucleotides within the sequences
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    complementDna = ''
    for index in range(len(dna)):
        complementDna = get_complement(dna[index]) + complementDna
    return complementDna


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        The added unit test is to test to make sure that the entire string is
        returned if there is no end condon
    >>> rest_of_ORF("ATGTGCCC")
    'ATGTGCCC'
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    sequence = 'ATG'
    index = 3
    while index < len(dna)-2:
        if dna[index:index+3] == 'TAG' or dna[index:index+3] == 'TAA' or dna[index:index+3] == 'TGA':
            return sequence
        sequence += dna[index:index+3]
        index += 3
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    the added test just provides anohter option to make sure that the function
    doesn't just work on this one specific example
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    >>> find_all_ORFs_oneframe("ATGCCCCCCCCCTGACGGATGCCCCGG")
    ['ATGCCCCCCCCC', 'ATGCCCCGG']
    """
    skipUntil = 0
    sequences = []
    temp = ''
    for index in range(len(dna)):
        if index < len(dna) and index >= skipUntil and index >= 3 and index % 3 == 0:
            if dna[index-3:index] == 'ATG':
                temp = rest_of_ORF(dna[index-3:])
                sequences.append(temp)
                skipUntil = index + len(temp)
    return sequences


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        the extra test was to make sure that the function works with more than
        just a single input
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']

    >>> find_all_ORFs("CATGCCATGCTGA")
    ['ATGCTGA', 'ATGCCATGC']

    """
    temp = []
    orfs = []
    for index in range(3):
        temp = find_all_ORFs_oneframe(dna[index:])
        for orf in temp:
            orfs.append(orf)
    return orfs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    strand1 = find_all_ORFs(dna)
    strand2 = find_all_ORFs(get_reverse_complement(dna))
    return strand1 + strand2


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
