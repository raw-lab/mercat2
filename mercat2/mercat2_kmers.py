# -*- coding: utf-8 -*-
"""mercat2_kmers.py: Module for calculating k-mer counts
"""

import humanize
import timeit
from collections import OrderedDict
from joblib import Parallel, delayed


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
def find_kmers(file:str, kmer:int, min_count:int, num_cores:int):
    '''Calculates the k-mer count in a fasta file.

    Parameters:
        file (str): path to a fasta file to scan for k-mers.
        kmer (int): k-mer length.
        min_count (int): minimum count of k-mers found to be considered significant.
        num_cores (int): cpu count to use for parallel processing.

    Returns:
        dict: A dictionary with the counts of each k-mer found.
    '''

    start_time = timeit.default_timer()

    print("Loading Sequences")
    sequences = OrderedDict()
    with open(file, 'r') as f:
        seq = ""
        sname = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sname:
                    sequences[sname] = ""
                if seq:
                    sequences[sname] = seq
                    seq = ""
                sname = line[1:]
                sname = sname.split("#",1)[0].strip()
            else:
                line = line.replace("*","")
                seq += line
        sequences[sname] = seq

    print(("Number of sequences in " + file + " = "+ str(humanize.intword(len(sequences)))))
    
    results = Parallel(n_jobs=num_cores)(
        delayed(calculateKmerCount)(sequences[seq], kmer) for seq in sequences)

    kmerlist = dict()
    for res in results:
        for k,v in res.items():
            if k in kmerlist:
                kmerlist[k] += v
            else: kmerlist[k] = v

    print(f"Time to compute {kmer}-mers: {round(timeit.default_timer() - start_time,2)} secs")

    significant_kmers = dict()
    for k,v in kmerlist.items():
        if kmerlist[k] >= min_count:
            significant_kmers[k] = v

    return significant_kmers
