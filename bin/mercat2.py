#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""mercat2.py: Python code for Parallel k-mer counting."""

__version__     = "0.3"
__author__      = "Jose L. Figueroa, Richard A. White III"
__copyright__   = "Copyright 2022"

import os
import psutil
from distutils.util import strtobool
import ray
import argparse
import timeit

# Mercat libraries
from mercat2_lib import (mercat2_fasta, mercat2_Chunker, mercat2_kmers, mercat2_diversity, mercat2_figures, mercat2_report)


# GLOBAL VARIABLES
FILE_EXT_FASTQ = ['.fq', '.fastq']
FILE_EXT_NUCLEOTIDE = [".fasta", ".fa", ".fna", ".ffn"]
FILE_EXT_PROTEIN = [".faa"]


# RAM Usage
def mem_use():
    return round(psutil.virtual_memory().used/1024.0**3, 2)


## Parse Arguments
#
def parseargs():
    '''Returns the parsed command line options.'''
    num_cores = psutil.cpu_count(logical=False)
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', type=str, required=False, help='path-to-input-file') #default=nucleotide
    parser.add_argument('-f', type=str, required=False, help='path-to-folder-containing-input-files')
    parser.add_argument('-k', type=int, required = True, help='kmer length')
    parser.add_argument('-n', type=int, default=num_cores, help='no of cores [auto detect]')  # no of cores to use
    parser.add_argument('-c', type=int, default=10, help='minimum kmer count [10]')  # minimum kmer count to report
    parser.add_argument('-prod', action='store_true', help='run Prodigal on fasta files')
    parser.add_argument('-fgs', action='store_true', help='run FragGeneScanRS on fasta files')
    parser.add_argument('-s', type=int, default=100, required=False, help='Split into x MB files. [100]')
    parser.add_argument('-o', type=str, default='mercat_results', required=False, help="Output folder, default = 'mercat_results' in current directory")
    parser.add_argument('-lowmem', type=strtobool, default=None, help="Flag to use incremental PCA when low memory is available. [auto]")
    parser.add_argument('-no_metaomestats', action='store_true', help='run prodigal on fasta files')
    parser.add_argument('-category_file', type=str, default=None, help=argparse.SUPPRESS)

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
        if not mercat2_fasta.check_command('prodigal'):
            exit(1)
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


## Create Figures
#
def createFigures(tsv_list:dict, type_string:str, out_path:os.PathLike, lowmem=None, class_file=None):
    print(f"\nCreating {type_string} Graphs")
    print(f"Virtual Memory {mem_use()}GB")

    figPlots = dict()
    start_time = timeit.default_timer()

    # Stacked Bar Plots (top kmer counts)
    combined_tsv = os.path.join(out_path, 'combined.tsv')
    if not os.path.exists(combined_tsv):
        mercat2_report.merge_tsv(tsv_list, combined_tsv)
        print(f"\nTime to merge TSV files: {round(timeit.default_timer() - start_time,2)} seconds")
        print(f"Virtual Memory {mem_use()}GB")
    
    figPlots[f"Combined {type_string} kmer Summary"] = mercat2_figures.kmer_summary(combined_tsv)
    print(f"Time to compute Combined {type_string} kmer Summary: {round(timeit.default_timer() - start_time,2)} seconds")
    print(f"Virtual Memory {mem_use()}GB")

    # PCA
    if len(tsv_list) > 3:
        print("\nRunning PCA")
        print(f"Virtual Memory {mem_use()}GB")
        start_time = timeit.default_timer()

        combined_tsv = os.path.join(out_path, 'combined_T.tsv')
        if not os.path.exists(combined_tsv):
            mercat2_report.merge_tsv_T(tsv_list, combined_tsv)
            print(f"\nTime to merge TSV files: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")

        figPlots[f"{type_string} PCA"] = mercat2_figures.plot_PCA(combined_tsv, out_path, lowmem, class_file)

        print(f"Virtual Memory {mem_use()}GB")

    return figPlots


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
    m_flag_fgs = __args__.fgs
    m_chunk_size = __args__.s
    m_outputfolder = __args__.o
    m_lowmem = None if __args__.lowmem is None else bool(__args__.lowmem)
    m_class_file = __args__.category_file
    
    #if os.path.exists(m_outputfolder) and os.path.isdir(m_outputfolder):
    #    shutil.rmtree(m_outputfolder)
    if not m_outputfolder:
        m_outputfolder = 'mercat_results'
    os.makedirs(m_outputfolder, exist_ok=True)

    print(f"\nVirtual Memory {mem_use()}GB")

    # Initialize Ray
    try:
        ray.init(address='auto', log_to_driver=False) # First try if ray is setup for a cluster
    except:
        ray.init(log_to_driver=False)

    # Load input files
    gc_content = dict()
    files_nucleotide = {}
    files_protein = {}
    cleanpath = os.path.join(m_outputfolder, 'clean')
    print("Loading files")
    start_time = timeit.default_timer()
    if m_inputfolder:
        m_inputfolder = os.path.abspath(os.path.expanduser(m_inputfolder))
        for fname in os.listdir(m_inputfolder):
            file = os.path.join(m_inputfolder, fname)
            if not os.path.isdir(file):
                basename, f_ext = os.path.splitext(fname)
                if f_ext in FILE_EXT_FASTQ:
                    files_nucleotide[basename] = mercat2_fasta.fastq_processing(file, cleanpath, basename)
                elif f_ext in FILE_EXT_NUCLEOTIDE:
                    file,stat = mercat2_fasta.removeN(file, cleanpath)
                    files_nucleotide[basename] = file
                    gc_content[basename] = stat['GC Content']
                    #TODO:MetaomeStats
                elif f_ext in FILE_EXT_PROTEIN:
                    files_protein[basename] = file
    else:
        basepath = os.path.abspath(os.path.expanduser(m_inputfile))
        basename,f_ext = os.path.splitext(os.path.basename(basepath))
        if f_ext in FILE_EXT_FASTQ:
            m_inputfile = mercat2_fasta.fastq_processing(m_inputfile, cleanpath, basename)
            files_nucleotide[basename] = m_inputfile
        elif f_ext in FILE_EXT_NUCLEOTIDE:
            m_inputfile,stat = mercat2_fasta.removeN(m_inputfile, cleanpath)
            files_nucleotide[basename] = m_inputfile
            gc_content[basename] = stat['GC Content']
            #TODO:MetaomeStats
        elif f_ext in FILE_EXT_PROTEIN:
            files_protein[basename] = m_inputfile
    print(f"Time to load {len(files_nucleotide)+len(files_protein)} files: {round(timeit.default_timer() - start_time,2)} seconds")

    # ORF Call if -prod flag
    if m_flag_prodigal:
        @ray.remote(num_cpus=1)
        def orf_call(basename, file, prodpath):
            return mercat2_fasta.orf_call(basename, file, prodpath)

        print(f"Running Prodigal on {len(files_nucleotide)} files")
        prodpath = os.path.join(m_outputfolder, 'prodigal')
        jobs = []
        for basename,file in files_nucleotide.items():
            jobs += [orf_call.remote(basename, file, prodpath)]
        while jobs:
            ready,jobs = ray.wait(jobs)
            name,amino = ray.get(ready[0])
            files_protein[name] = amino

    # ORF Call if -fgs flag
    if m_flag_fgs:
        @ray.remote(num_cpus=1)
        def orf_call_fgs(basename, file, prodpath):
            return mercat2_fasta.orf_call_fgs(basename, file, outpath)

        print(f"Running FragGeneScanRS on {len(files_nucleotide)} files")
        outpath = os.path.join(m_outputfolder, 'fgs')
        jobs = []
        for basename,file in files_nucleotide.items():
            jobs += [orf_call_fgs.remote(basename, file, outpath)]
        while jobs:
            ready,jobs = ray.wait(jobs)
            ret = ray.get(ready[0])
            if ret:
                files_protein[ret[0]] = ret[1]

    # Chunk large files
    if m_chunk_size > 0:
        print("Checking for large files")
        start_time = timeit.default_timer()
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
    print(f"Time to check for large files: {round(timeit.default_timer() - start_time,2)} seconds")

    # Begin processing files
    figPlots = dict()
    tsv_stats = dict()
    os.makedirs(os.path.join(m_outputfolder, 'tsv'), exist_ok=True)

    @ray.remote(num_cpus=1)
    def countKmers(file):
        return mercat2_kmers.find_kmers(file, m_kmer, m_min_count, m_num_cores)

    @ray.remote(num_cpus=1)
    def run_mercat2(basename:str, files:list, out_file:os.PathLike):
        #start_time = timeit.default_timer()
        kmers = dict()
        jobs = []
        for file in files:
            jobs += [countKmers.remote(file)]
        while(jobs):
            ready,jobs = ray.wait(jobs)
            for k,v in ray.get(ready[0]).items():
                if k in kmers:
                    kmers[k] += v
                else:
                    kmers[k] = v
        #print(f"Time to compute {m_kmer}-mers: {round(timeit.default_timer() - start_time,2)} secs")
        if len(kmers):
            print(f"Significant k-mers: {len(kmers)}")
            with open(out_file, 'w') as writer:
                print('k-mer', f'{basename}_Count', sep='\t', file=writer)
                for kmer,count in sorted(kmers.items()):
                    print(kmer, count, sep='\t', file=writer)
            return basename,out_file
        else:
            print(f"No significant k-mers found")
            return None


    ## Process Nucleotides ##
    if len(files_nucleotide):
        print("Processing Nucleotides")
        print(f"Running Mercat2 using {m_num_cores} cores")
    jobs = []
    start_time = timeit.default_timer()
    for basename,files in files_nucleotide.items():
        out_counts = out_counts = os.path.join(m_outputfolder, 'tsv', f"{basename}_nucleotide_counts.tsv")
        if not os.path.exists(out_counts):
            jobs += [run_mercat2.remote(basename, files, out_counts)]
    df_list = dict()
    while(jobs):
        ready,jobs = ray.wait(jobs)
        if ready[0]:
            try:
                basename,kmers = ray.get(ready[0])
                df_list[basename] = kmers
            except: pass
    if len(files_nucleotide):
        print(f"Time to compute {m_kmer}-mers: {round(timeit.default_timer() - start_time,2)} seconds")
    # Stacked Bar Plots (top kmer counts)
    if len(df_list):
        figPlots.update(createFigures(df_list, "Nucleotide", m_outputfolder, m_lowmem, m_class_file))
    if len(gc_content) > 0:
        figPlots['Sample GC Summary'] = mercat2_figures.GC_plot_sample(gc_content)


    ## Process Proteins ##
    if len(files_protein):
        print("Processing Proteins")
        print(f"Running Mercat2 using {m_num_cores} cores")
    jobs = []
    start_time = timeit.default_timer()
    df_list = dict()
    for basename,files in files_protein.items():
        out_counts = out_counts = os.path.join(m_outputfolder, 'tsv', f"{basename}_protein_counts.tsv")
        if os.path.exists(out_counts):
            df_list[basename] = out_counts
        else:
            jobs += [run_mercat2.remote(basename, files, out_counts)]
    while(jobs):
        ready,jobs = ray.wait(jobs)
        if ready[0]:
            try:
                basename,kmers = ray.get(ready[0])
                df_list[basename] = kmers
            except: pass
    if len(files_protein):
        print(f"Time to compute {m_kmer}-mers: {round(timeit.default_timer() - start_time,2)} seconds")
    if len(df_list):
        figPlots.update(createFigures(df_list, "Protein", m_outputfolder, m_lowmem, m_class_file))

    # Plot Data
    report_dir = os.path.join(m_outputfolder, "report")
    os.makedirs(report_dir, exist_ok=True)

    @ray.remote(num_cpus=1)
    def diversity(infile, outfile):
        mercat2_diversity.compute_alpha_beta_diversity(infile, outfile)
        return
    jobs = list()
    for filename in os.listdir(os.path.join(m_outputfolder, 'tsv')):
        if not filename.endswith('.tsv'):
            continue
        jobs += [diversity.remote(os.path.join(m_outputfolder, 'tsv', filename), os.path.join(report_dir, f'diversity-{filename}'))]

    while jobs:
        ready,jobs = ray.wait(jobs)

    mercat2_report.write_html(os.path.join(report_dir, "report.html"), figPlots, tsv_stats)
    figPlots = mercat2_figures.plot_sample_metrics(files_protein, report_dir)
    mercat2_report.write_html(os.path.join(report_dir, "protein_metrics.html"), figPlots, dict())
    return


# If main app
if __name__ == "__main__":
    mercat_main()
