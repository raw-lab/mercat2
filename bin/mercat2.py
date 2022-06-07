#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""mercat2.py: Python code for Parallel k-mer counting."""

__version__     = "0.2"
__author__      = "Jose L. Figueroa III, Richard A. White III"
__copyright__   = "Copyright 2022"

import sys
import os
import psutil
import shutil
import subprocess
import ray
import pandas as pd
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import timeit

# Mercat libraries
from mercat2_lib import (mercat2_fasta, mercat2_Chunker, mercat2_kmers, mercat2_report)


# GLOBAL VARIABLES
FILE_EXT_FASTQ = ['.fastq']
FILE_EXT_NUCLEOTIDE = [".fasta", ".fa", ".fna", ".ffn"]
FILE_EXT_PROTEIN = [".faa"]


## Check Command
#
def check_command(cmd:str):
    '''Raises an error and exists if the given command fails'''
    try:
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Mercat Error: {cmd} not found, please setup {cmd} using: 'conda install {cmd}'")
        sys.exit(1)
    return


## Parse Arguments
#
def parseargs():
    '''Returns the parsed command line options.'''
    num_cores = psutil.cpu_count(logical=False)
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-i', type=str, required=False, help='path-to-input-file') #default=nucleotide
    parser.add_argument('-f', type=str, required=False, help='path-to-folder-containing-input-files')
    parser.add_argument('-k', type=int, required = True, help='kmer length')
    parser.add_argument('-n', type=int, default=num_cores, help='no of cores [default = all]')  # no of cores to use
    parser.add_argument('-c', type=int, default=10, help='minimum kmer count [default = 10]')  # minimum kmer count to report
    parser.add_argument('-prod', action='store_true', help='run prodigal on fasta file')
    parser.add_argument('-s', type=int, default=100, required=False, help='Split into x MB files. Default = 100MB')
    parser.add_argument('-o', type=str, default='mercat_results', required=False, help="Output folder, default = 'mercat_results' in current directory")

    # Process arguments
    args = parser.parse_args()

    # Check folder/input file options
    if args.i and args.f:
        parser.error("Can only specify either an input file (-i) or path to folder containing input files (-f) at a time")
    if args.i:
        if not os.path.isfile(args.i):
            parser.error(f"file {args.i} is not valid.\n")
    elif args.f:
        if not os.path.isdir(args.f):
            parser.error(f"folder {args.f} is not valid.\n")
    else:
        parser.error("Please provide either an input file (-i) or an input folder (-f)")

    # check prodigal/protein flags
    if args.prod:
        check_command('prodigal')
    return args


## Chunk Files
#
@ray.remote(num_cpus=1)
def chunk_files(name:str, file:str, chunk_size:int, outpath:str):
    '''Checks if the given file is larger than chunk_size and splits the file as necessary.
    Fasta files are split at sequence headers to keep individual sequences contiguous.

    Parameters:
        name (str): The name of the sample file for proper grouping of the files.
        file (str): Path to a fasta.
        chunk_size (int): The size, in MB, to split the files at.
        outpath (str): Path to a folder to save the chunked files in. Folder should not exist.

    Returns:
        tuple: A tuple with the name, and list of paths to the chunked files.
    '''

    if os.stat(file).st_size >= (chunk_size*1024*1024):
        print(f"Large input file provided: {file}\nSplitting it into smaller files...\n")
        os.makedirs(outpath, exist_ok=True)
        all_chunks = mercat2_Chunker.Chunker(file, outpath, str(chunk_size)+"M", ">").files
    else:
        all_chunks = [file]
    return (name, all_chunks)


## Main
#
def mercat_main():
    '''Main method'''

    # Process Arguments
    __args__ = parseargs()
    # Set Flags
    m_kmer = __args__.k
    m_num_cores = __args__.n
    m_inputfile = __args__.i
    m_inputfolder = __args__.f
    m_min_count = __args__.c
    m_flag_prodigal = __args__.prod
    m_chunk_size = __args__.s
    m_outputfolder = __args__.o
    
    if os.path.exists(m_outputfolder) and os.path.isdir(m_outputfolder):
        shutil.rmtree(m_outputfolder)
    os.makedirs(m_outputfolder, exist_ok=True)

    # Initialize Ray
    try:
        ray.init(address='auto', log_to_driver=False) # First try if ray is setup for a cluster
    except:
        ray.init(log_to_driver=False)

    # Load input files
    all_ipfiles = []
    gc_content = dict()
    files_nucleotide = {}
    files_protein = {}
    cleanpath = os.path.join(m_outputfolder, 'clean')
    print("Loading files")
    if m_inputfolder:
        m_inputfolder = os.path.abspath(os.path.expanduser(m_inputfolder))
        for fname in os.listdir(m_inputfolder):
            file = os.path.join(m_inputfolder, fname)
            if not os.path.isdir(file):
                # skip directories
                basename, f_ext = os.path.splitext(fname)
                if f_ext in FILE_EXT_FASTQ:
                    check_command('fastqc')
                    check_command('fastp')
                    file = mercat2_fasta.fastq_processing(file, cleanpath, basename)
                    files_nucleotide[basename] = file
                elif f_ext in FILE_EXT_NUCLEOTIDE:
                    file,stat = mercat2_fasta.removeN(file, cleanpath)
                    files_nucleotide[basename] = file
                    gc_content[basename] = stat['GC Content']
                elif f_ext in FILE_EXT_PROTEIN:
                    files_protein[basename] = file
                all_ipfiles.append(file)
    else:
        basepath = os.path.abspath(os.path.expanduser(m_inputfile))
        basename, f_ext = os.path.splitext(os.path.basename(basepath))
        if f_ext in FILE_EXT_FASTQ:
            check_command('fastqc')
            check_command('fastp')
            m_inputfile = mercat2_fasta.fastq_processing(m_inputfile, cleanpath, basename)
            files_nucleotide[basename] = m_inputfile
        elif f_ext in FILE_EXT_NUCLEOTIDE:
            m_inputfile,stat = mercat2_fasta.removeN(m_inputfile, cleanpath)
            files_nucleotide[basename] = m_inputfile
            gc_content[basename] = stat['GC Content']
        elif f_ext in FILE_EXT_PROTEIN:
            files_protein[basename] = m_inputfile
        all_ipfiles.append(m_inputfile)
        all_ipfiles.append(os.path.abspath(m_inputfile))

    
    # ORF Call if -prod flag
    if m_flag_prodigal:
        @ray.remote(num_cpus=1)
        def orf_call(basename, file, prodpath):
            return mercat2_fasta.orf_call(basename, file, prodpath)

        print(f"Running prodigal on {len(files_nucleotide)} files")
        prodpath = os.path.join(m_outputfolder, 'prodigal')
        jobs = []
        for basename,file in files_nucleotide.items():
            jobs += [orf_call.remote(basename, file, prodpath)]
        while jobs:
            ready,jobs = ray.wait(jobs)
            name,amino = ray.get(ready[0])
            files_protein[name] = amino


    # Chunk large files
    if m_chunk_size > 0:
        print("Checking for large files")
        dir_chunks = os.path.join(m_outputfolder, "chunks")
        # nucleotides
        jobs = []
        for basename,file in files_nucleotide.items():
            chunk_path = os.path.join(dir_chunks, f"{basename}_nucleotide")
            jobs += [chunk_files.remote(basename, file, m_chunk_size, chunk_path)]
        while jobs:
            ready,jobs = ray.wait(jobs)
            name,chunks = ray.get(ready[0])
            files_nucleotide[name] = chunks
        # proteins
        jobs = []
        for basename,file in files_protein.items():
            chunk_path = os.path.join(dir_chunks, f"{basename}_protein")
            jobs += [chunk_files.remote(basename, file, m_chunk_size, chunk_path)]
        while jobs:
            ready,jobs = ray.wait(jobs)
            name,chunks = ray.get(ready[0])
            files_protein[name] = chunks

    # Begin processing files
    figPlots = dict()
    tsv_stats = dict()
    os.makedirs(os.path.join(m_outputfolder, 'tsv'), exist_ok=True)

    @ray.remote(num_cpus=1)
    def countKmers(file, kmer, min_count, num_cores):
        return mercat2_kmers.find_kmers(file, kmer, min_count, num_cores)

    # Process Nucleotides
    if len(files_nucleotide):
        print("Processing Nucleotides")
    df_list = dict()
    for basename,files in files_nucleotide.items():
        np_string = '_nucleotide'

        # Run Mercat
        print(f"Running mercat using {m_num_cores} cores on {basename}{np_string}")
        start_time = timeit.default_timer()
        kmers = dict()
        jobs = []
        for file in files:
            jobs += [countKmers.remote(file, m_kmer, m_min_count, m_num_cores)]
        while(jobs):
            ready,jobs = ray.wait(jobs)
            for k,v in ray.get(ready[0]).items():
                if k in kmers:
                    kmers[k] += v
                else:
                    kmers[k] = v
        print(f"Time to compute {m_kmer}-mers: {round(timeit.default_timer() - start_time,2)} secs")
        if len(kmers):
            print(f"Significant k-mers: {len(kmers)}")
            df_list[basename] = pd.DataFrame(index=kmers.keys(), data=kmers.values(), columns=['Count'])
            outfile = os.path.join(m_outputfolder, 'tsv', f"{basename}{np_string}_summary.tsv")
            dfGC = df_list[basename].reset_index().rename(columns=dict(index='k-mer'))
            dfGC['GC%'] = dfGC['k-mer'].apply(lambda k: round(100.0 * (k.count('G') + k.count('C')) / len(k), 2))
            dfGC.to_csv(outfile, index=False, sep='\t')
        else:
            print(f"No significant k-mers found")
    # Stacked Bar Plots (top kmer counts)
    if len(df_list):
        figPlots["Combined Nucleotide kmer Summary"] = mercat2_report.kmer_summary(df_list)
        if len(files_nucleotide) > 3:
            dfPCA = pd.DataFrame()
            for name,df in df_list.items():
                dfTemp = df.rename(columns=dict(Count=name))
                dfPCA = dfPCA.merge(dfTemp, left_index=True, right_index=True, how="outer")
            dfPCA = dfPCA.reset_index().rename(columns=dict(index="k-mer"))
            figPlots["Nucleotide PCA"] = mercat2_report.PCA_plot(dfPCA)
        if m_kmer >= 10:
            figPlots['k-mer GC Summary'] = mercat2_report.GC_plot_kmer(df_list)
    if len(gc_content) > 0:
        figPlots['Sample GC Summary'] = mercat2_report.GC_plot_sample(gc_content)
    

    # Process Proteins
    if len(files_protein):
        print("Processing Proteins")
    df_list = dict()
    for basename,files in files_protein.items():
        np_string = '_protein'

        # Run Mercat
        start_time = timeit.default_timer()
        kmers = dict()
        jobs = []
        for file in files:
            jobs += [countKmers.remote(file, m_kmer, m_min_count, m_num_cores)]
        while(jobs):
            ready,jobs = ray.wait(jobs)
            for k,v in ray.get(ready[0]).items():
                if k in kmers:
                    kmers[k] += v
                else:
                    kmers[k] = v
        print(f"Time to compute {m_kmer}-mers: {round(timeit.default_timer() - start_time,2)} secs")
        if len(kmers):
            print("Significant k-mers:", len(kmers))
            df_list[basename] = pd.DataFrame(index=kmers.keys(), data=kmers.values(), columns=['Count'])
            outfile = os.path.join(m_outputfolder, 'tsv', f"{basename}{np_string}_summary.tsv")
            df_list[basename].to_csv(outfile, index=False, sep='\t')
        else:
            print(f"No significant k-mers found")
    # Stacked Bar Plots (top kmer counts)
    if len(df_list):
        figPlots["Combined Protein kmer Summary"] = mercat2_report.kmer_summary(df_list)
        if len(files_protein) > 3:
            dfPCA = pd.DataFrame()
            for name,df in df_list.items():
                df = df.rename(columns=dict(Count=name))
                dfPCA = dfPCA.merge(df, left_index=True, right_index=True, how="outer")
            dfPCA = dfPCA.reset_index().rename(columns=dict(index="k-mer"))
            figPlots["Protein PCA"] = mercat2_report.PCA_plot(dfPCA)
    if len(files_protein):
        figPlots["Sample Protein Metrics Summary"],tsv_stats["Protein_Metrics"] = mercat2_report.plot_sample_metrics(files_protein)
    
    
    # Plot Data
    plots_dir = os.path.join(m_outputfolder, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    mercat2_report.write_html(os.path.join(plots_dir, "report.html"), figPlots, tsv_stats)
    return


# If main app
if __name__ == "__main__":
    mercat_main()
