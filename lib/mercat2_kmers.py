# -*- coding: utf-8 -*-
"""mercat2_kmers.py: Module for calculating k-mer counts
"""

import gzip
from pathlib import Path

## Calculate k-mer count
#
def calculateKmerCount(seq:str, kmer:int) -> dict:
    '''Calculates the k-mer count in a sequence.
    This is a helper method for find_kmers()

    Parameters:
        seq (str): The sequence to search for k-mers.
        kmer (int): The k-mer size

    Returns:
        dict: A dictionary with the counts of each k-mer found.
    '''

    kmerlist = dict()
    for i in range((len(seq) - kmer) + 1):
        k = seq[i:i+kmer]
        if k not in kmerlist:
            kmerlist[k] = 0
        kmerlist[k] += 1
    return kmerlist


## Find k-mers
def find_kmers(file:Path, kmer:int, min_count:int):
    '''Calculates the k-mer count in a fasta file.

    Parameters:
        file (Path): path to a fasta file to scan for k-mers.
        kmer (int): k-mer length.
        min_count (int): minimum count of k-mers found to be considered significant.

    Returns:
        dict: A dictionary with the counts of each k-mer found.
    '''

    num = 0
    kmerlist = dict()

    f = gzip.open(file, 'rt') if file.suffix=='.gz' else open(file, 'r')

    seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            num += 1
            if seq:
                # count k-mers
                for i in range((len(seq) - kmer) + 1):
                    k = seq[i:i+kmer]
                    if k not in kmerlist:
                        kmerlist[k] = 0
                    kmerlist[k] += 1
                seq = ""
        else:
            seq += line.replace("*","")
    # Count last sequence
    for i in range((len(seq) - kmer) + 1):
        k = seq[i:i+kmer]
        if k not in kmerlist:
            kmerlist[k] = 0
        kmerlist[k] += 1

    f.close()

    significant_kmers = dict()
    for k,v in kmerlist.items():
        if kmerlist[k] >= min_count:
            significant_kmers[k] = v

    return significant_kmers

