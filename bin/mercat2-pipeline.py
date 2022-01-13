#!/usr/bin/env python

"""mercat2.py: Python code for Parallel k-mer counting."""

__author__      = "Jose L. Figueroa III, Richard A. White III"
__copyright__   = "Copyright 2022"

import sys
import os
import psutil
import shutil
import subprocess
import pandas as pd
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from joblib import Parallel, delayed

import dask.dataframe as dd

# Mercat libraries
from mercat2 import (mercat2_fasta, mercat2_Chunker, mercat2_kmers, mercat2_report, mercat2_metrics)


# GLOBAL VARIABLES
FILE_EXT_FASTQ = ['.fastq', '.fastq.gz']
FILE_EXT_NUCLEOTIDE = [".fasta", ".fa", ".fna", ".ffn", ".fasta.gz", ".fa.gz", ".fna.gz", ".ffn.gz"]
FILE_EXT_PROTEIN = [".faa", ".faa.gz"]


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
    parser.add_argument('-prot', action='store_true', help='assume fasta files are protein input files')
    parser.add_argument('-prod', action='store_true', help='run prodigal on fasta file')
    parser.add_argument('-frag', action='store_true', help='run FragGeneScan+ on fasta file')
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
        if args.prot:
            parser.error("Can only provide one of -prot or -prod option at a time")
        check_command('prodigal')
    return args
def check_args(ipfile,args,def_option,m_parser): #TODO: Merge into parseargs()
    given_ext = (os.path.splitext(ipfile)[1]).strip()
    if def_option:
        if given_ext not in FILE_EXT_PROTEIN:
            m_parser.error("Input file provided should be one of the following formats: " + str(FILE_EXT_PROTEIN))

    if args.prod:
        # if not args.q and given_ext not in protein_file_ext:
        if given_ext not in FILE_EXT_PROTEIN:
            m_parser.error("Input file provided should be one of the following formats: " + str(FILE_EXT_PROTEIN))
    
    if args.prot:
        if given_ext != ".faa":
            m_parser.error("Input file provided should be in .faa format")


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
    m_flag_protein = __args__.prot
    m_chunk_size = __args__.s
    m_kmerstring = f"{m_kmer}-mers"
    m_outputfolder = __args__.o
    os.makedirs(m_outputfolder, exist_ok=True)

    np_string = 'protein' if m_flag_protein or m_flag_prodigal else 'nucleotide' #TODO:Rename variable

    # Load input files
    all_ipfiles = []
    gc_content = dict()
    files_nucleotide = {}
    files_protein = {}
    cleanpath = os.path.join(m_outputfolder, 'clean')
    print("Loading files")
    if m_inputfolder:
        m_inputfolder = os.path.abspath(m_inputfolder)
        for fname in os.listdir(m_inputfolder):
            file = os.path.join(m_inputfolder, fname)
            if not os.path.isdir(file):
                # skip directories
                basename, f_ext = os.path.splitext(fname)
                if f_ext in FILE_EXT_FASTQ:
                    print("***FASTQ PROCESS***")
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
        basename, f_ext = os.path.splitext(m_inputfile)
        if f_ext in FILE_EXT_FASTQ:
            print("***FASTQ PROCESS***")
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
        print(f"Running prodigal on {len(files_nucleotide)} files")
        prodpath = os.path.join(m_outputfolder, 'prodigal')
        results = Parallel(n_jobs=m_num_cores)(
            delayed(mercat2_fasta.orf_call)(basename, file, prodpath) for basename,file in files_nucleotide.items())
        for res in results:
            for name,file in res.items():
                files_protein[name] = file

    # Begin processing files
    figPlots = dict()
    tsv_stats = dict()
    os.makedirs(os.path.join(m_outputfolder, 'tsv'), exist_ok=True)

    # Process Nucleotides
    print("Processing Nucleotides")
    df_list = dict()
    for basename, file in files_nucleotide.items():
        print(basename, file)
        np_string = '_nucleotide'

        # Split into chunks
        file_size = os.stat(file).st_size
        if file_size >= (m_chunk_size*1024*1024):
            print("Large input file provided: Splitting it into smaller files...\n")
            dir_chunks = os.path.join(m_outputfolder, "chunks", f"{basename}_{np_string}")
            os.makedirs(dir_chunks, exist_ok=True)
            all_chunks = mercat2_Chunker.Chunker(file, dir_chunks, str(m_chunk_size)+"M", ">").files
        else:
            all_chunks = [file]

        # Run Mercat
        print(("Running mercat using " + str(m_num_cores) + " cores"))
        print(("input file: " + file))
        for chunk in all_chunks:
            df = mercat2_kmers.find_kmers(chunk, m_kmer, m_min_count, m_num_cores)
            df_list[basename] = df
            outfile = os.path.join(m_outputfolder, 'tsv', f"{basename}{np_string}_summary.tsv")
            df.to_csv(outfile, index_label="k-mers", index=True, sep='\t')
    # Stacked Bar Plots (top kmer counts)
    if len(df_list):
        figPlots["Combined Nucleotide kmer Summary"] = mercat2_report.kmer_summary(df_list)
        if len(files_nucleotide) > 3:
            dfPCA = pd.DataFrame()
            for name,df in df_list.items():
                dfTemp = df.rename(columns=dict(Count=name))
                dfPCA = dfPCA.merge(dfTemp, left_index=True, right_index=True, how="outer")
            dfPCA = dfPCA.reset_index().rename(columns=dict(index="k-mer"))
            figPlots["Nucleotide PCA"] = [mercat2_report.PCA_plot(dfPCA)]
        if m_kmer >= 10:
            figPlots['k-mer GC Summary'] = [mercat2_report.GC_plot_kmer(df_list)]
    if len(gc_content) > 0:
        figPlots['Sample GC Summary'] = [mercat2_report.GC_plot_sample(gc_content)]
    

    # Process Proteins
    print("Processing Proteins")
    df_prot_stats = pd.DataFrame(columns=['sample', 'PI', 'MW', 'Hydro'])
    df_list = dict()
    for basename, file in files_protein.items():
        print(basename, file)
        np_string = '_protein'

        if file_size >= (m_chunk_size*1024*1024):
            print("Large input file provided: Splitting it into smaller files...\n")
            dir_chunks = os.path.join(m_outputfolder, "chunks", f"{basename}_{np_string}")
            os.makedirs(dir_chunks, exist_ok=True)
            all_chunks = mercat2_Chunker.Chunker(file, dir_chunks, str(m_chunk_size)+"M", ">").files
        else:
            all_chunks = [file]

        # Run Mercat
        print(("Running mercat using " + str(m_num_cores) + " cores"))
        print(("input file: " + file))
        for chunk in all_chunks:
            df = mercat2_kmers.find_kmers(chunk, m_kmer, m_min_count, m_num_cores)
            df_list[basename] = df
            outfile = os.path.join(m_outputfolder, 'tsv', f"{basename}{np_string}_summary.tsv")
            df.to_csv(outfile, index_label="k-mers", index=True, sep='\t')
    # Stacked Bar Plots (top kmer counts)
    if len(df_list):
        figPlots["Combined Protein kmer Summary"] = mercat2_report.kmer_summary(df_list)
        if len(files_protein) > 3:
            dfPCA = pd.DataFrame()
            for name,df in df_list.items():
                df = df.rename(columns=dict(Count=name))
                dfPCA = dfPCA.merge(df, left_index=True, right_index=True, how="outer")
            dfPCA = dfPCA.reset_index().rename(columns=dict(index="k-mer"))
            figPlots["Protein PCA"] = [mercat2_report.PCA_plot(dfPCA)]


    # Plot Data
    plots_dir = os.path.join(m_outputfolder, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    mercat2_report.write_html(os.path.join(plots_dir, "report.html"), figPlots, tsv_stats)
    return


    # PCA
    subdir = m_inputfolder+"/mercat_results/"
    if len(all_ipfiles) >= 4:
       figPlots['PCA'] = mercat2_report.PCA_plot(subdir)

    # Metastat Plots
    if len(gc_content) > 0:
        figPlots['Sample GC Summary'] = mercat2_report.GC_plot_sample(gc_content)

    if m_flag_protein:
        pass
        #figPlots[basename_ipfile+'_PI'] = mercat2_report.scatter_plots(basename_ipfile, 'PI', top10_all_samples, kmerstring)
        #figPlots[basename_ipfile+'_MW'] = mercat2_report.scatter_plots(basename_ipfile, 'MW', df_top10, kmerstring)
        #figPlots[basename_ipfile+'_Hydro'] = mercat2_report.scatter_plots(basename_ipfile, 'Hydro', df_top10, kmerstring)
    else:
        if m_kmer >= 10:
            figPlots['k-mer GC Summary'] = mercat2_report.GC_plot_kmer(df_list)
        #figPlots[basename_ipfile+'_GC'] = mercat2_report.scatter_plots(basename_ipfile, 'GC_Percent', df_top10, kmerstring)
        #figPlots[basename_ipfile+'_AT'] = mercat2_report.scatter_plots(basename_ipfile, 'AT_Percent', df_top10, kmerstring)

    # Save HTML Report
    mercat2_report.write_html(os.path.join(plots_dir, "report.html"), figPlots, tsv_stats)
    return 0


if __name__ == "__main__":
    mercat_main()
