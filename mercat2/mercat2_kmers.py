# -*- coding: utf-8 -*-

import os
import shutil
import humanize
import subprocess
import timeit
from numpy import number
import pandas as pd
from collections import OrderedDict
from joblib import Parallel, delayed

from mercat2 import (mercat2_Chunker, mercat2_metrics)

## Calculate k-mer count
#
def calculateKmerCount(seq:str, kmer:int) -> dict:
    '''Calculates the k-mer count in a sequence.
    This is a helper method for find_kmers()

    Parameters:
        seq (str): The sequence to search for kmers.
        kmer (int): The k-mer size

    Returns:
        dict: A dictionary with the counts of each kmer found.
    '''

    kmerlist = dict()
    for i in range((len(seq) - kmer) + 1):
        k = seq[i:i+kmer]
        if k not in kmerlist:
            kmerlist[k] = 0
        kmerlist[k] += 1
    return kmerlist


## Find k-mers
#
def find_kmers(file:str, kmer:int, min_count, num_cores:int):
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

    df = pd.DataFrame(index=significant_kmers.keys(), data=significant_kmers.values(), columns=['Count'])
    return df
    significant_kmers = []
    for k in kmerlist:
        if kmerlist[k] >= min_count:
            significant_kmers.append(k)

    print(f"Total number of {kmer}-mers found: {humanize.intword(len(kmerlist))}")
    #print(m_kmerstring +  " with count >= " + str(m_top_count) + ": " + str(humanize.intword(len(significant_kmers))))

    if is_protein:
        df = pd.DataFrame(0.0, index=significant_kmers, columns=['Count',"PI","MW","Hydro"])
        for k in significant_kmers:
            df.at[k,'Count'] = kmerlist[k]
            df.at[k,'PI'] = mercat2_metrics.predict_isoelectric_point_ProMoST(k)
            df.at[k,'MW'] = mercat2_metrics.calculate_MW(k)
            df.at[k,'Hydro'] = mercat2_metrics.calculate_hydro(k)
    else:
        df = pd.DataFrame(0, index=significant_kmers, columns=['Count',"GC_Percent","AT_Percent"])
        for k in significant_kmers:
            c_kmer = k
            df.at[k,'Count'] = kmerlist[k]
            len_cseq = float(len(c_kmer))
            df.at[k,'GC_Percent'] = round(((c_kmer.count("G")+c_kmer.count("C")) / len_cseq) * 100.0)
            df.at[k,'AT_Percent'] = round(((c_kmer.count("A")+c_kmer.count("T")) / len_cseq) * 100.0)

    df.to_csv(basename_ipfile + "_summary.csv", index_label="k-mers", index=True)

    splitSummaryFiles = []
    splitSummaryFiles.append(basename_ipfile + "_summary.csv")
    print(f"Total time: {round(timeit.default_timer() - start_time,2)} secs")

    return df
    num_chunks = len(all_chunks_ipfile)
    df = dd.read_csv(splitSummaryFiles)
    dfgb = df.groupby("k-mers").sum()
    try:
        df_top10 = dfgb.nlargest(10, 'Count').compute()
    except:
        return None
    dfsum = dfgb.sum(0).compute()

    dfgb.to_csv("./" + basename_ipfile + "_finalSummary*.csv", index_label=m_kmerstring, name_function=lambda name: str(name))

    if m_flag_protein:
        df_top10[['PI', 'MW', 'Hydro']] = df_top10[['PI', 'MW', 'Hydro']] / num_chunks
    else:
        df_top10[['GC_Percent', 'AT_Percent']] = df_top10[['GC_Percent', 'AT_Percent']] / num_chunks

    top10_all_samples[sample_name] = [df_top10, dfsum.Count]

    all_counts = dfgb.Count.values.compute().astype(int)
    mercat2_metrics.compute_alpha_beta_diversity(all_counts, basename_ipfile)

    if is_chunked:
        for tempfile in all_chunks_ipfile:
            os.remove(tempfile)
        for sf in splitSummaryFiles:
            os.remove(sf)
    return df
