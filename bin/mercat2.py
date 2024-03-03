#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""versatile k-mer counter and diversity estimator for database independent property analysis (DIPA) for multi-omic analysis"""

__version__     = "1.4.0"
__author__      = "Jose L. Figueroa III, Richard A. White III"
__copyright__   = "Copyright (c) 2022-2024"

import os
import re
from pathlib import Path
import subprocess
import shutil
import psutil
from distutils.util import strtobool
import ray
import argparse
import timeit

# Mercat libraries
from mercat2_lib import (mercat2_fasta, mercat2_Chunker, mercat2_kmers, mercat2_diversity, mercat2_figures, mercat2_report)


# GLOBAL VARIABLES
FILE_EXT_FASTQ = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]
FILE_EXT_NUCLEOTIDE = [".fasta", ".fa", ".fna", ".ffn", ".fasta.gz", ".fa.gz", ".fna.gz", ".ffn.gz"]
FILE_EXT_PROTEIN = [".faa", ".faa.gz"]

# RAM Usage
def mem_use():
    return round(psutil.virtual_memory().used/1024.0**3, 2)


## Parse Arguments
#
def parseargs():
    '''Returns the parsed command line options.'''
    num_cores = psutil.cpu_count(logical=False)
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', required=False, default=list(), help='path to input file', nargs='+')
    parser.add_argument('-f', type=str, required=False, help='path to folder containing input files')
    parser.add_argument('-k', type=int, required = True, help='kmer length')
    parser.add_argument('-n', type=int, default=num_cores, help='no of cores [auto detect]')
    parser.add_argument('-c', type=int, default=10, help='minimum kmer count [10]')
    parser.add_argument('-prod', action='store_true', help='run Prodigal on fasta files')
    parser.add_argument('-fgs', action='store_true', help='run FragGeneScanRS on fasta files')
    parser.add_argument('-s', type=int, default=100, required=False, help='Split into x MB files. [100]')
    parser.add_argument('-o', type=str, default='mercat_results', required=False, help="Output folder, default = 'mercat_results' in current directory")
    parser.add_argument('-replace', action='store_true', help="Replace existing output directory [False]")
    parser.add_argument('-lowmem', type=strtobool, default=None, help="Flag to use incremental PCA when low memory is available. [auto]")
    parser.add_argument('-skipclean', action='store_true', help='skip trimming of fastq files')
    parser.add_argument('-toupper', action='store_true', help='convert all input sequences to uppercase')
    parser.add_argument('-category_file', type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument('-debug', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('-pca', action='store_true', help='create interactive PCA plot of the samples (minimum of 4 fasta files required)')
    parser.add_argument('--version', '-v', action='version',
                        version=f'MerCat2: \n version: {__version__}',
                        help='show the version number and exit')

    # Process arguments
    args = parser.parse_args()

    global DEBUG
    DEBUG = args.debug

    # Check folder/input file options
    if not args.i and not args.f:
        parser.error("Please provide either an input file (-i) or an input folder (-f)")
    for filename in args.i:
        if not os.path.isfile(filename):
            parser.error(f"file '{filename}' is not valid.\n")
    if args.f and not os.path.isdir(args.f):
        parser.error(f"folder {args.f} is not valid.\n")

    # check prodigal/protein flags
    if args.prod:
        if not mercat2_fasta.check_command('prodigal'):
            exit(1)

    return args, parser


## Chunk Files
#
@ray.remote(num_cpus=1)
def chunk_files(name:str, filename:str, chunk_size:int, outpath:str):
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

    if os.stat(filename).st_size >= (chunk_size*1024*1024):
        os.makedirs(outpath, exist_ok=True)
        all_chunks = mercat2_Chunker.Chunker(filename, outpath, str(chunk_size)+"M", ">").files
    else:
        all_chunks = [filename]
    return (name, all_chunks)

@ray.remote(num_cpus=1)
def diversity(key, infile, outfile, sample_type):
    mercat2_diversity.compute_alpha_diversity(key, infile, outfile)
    return (sample_type, outfile)
@ray.remote(num_cpus=1)
def countKmers(file, kmer, min_count):
    return mercat2_kmers.find_kmers(Path(file), kmer, min_count)
@ray.remote(num_cpus=1)
def run_mercat2(basename:str, files:list, out_file:os.PathLike, kmer, min_count, num_cores):
    kmers = dict()
    jobs = []
    for file in files:
        jobs += [countKmers.remote(file, kmer, min_count)]
    while(jobs):
        ready,jobs = ray.wait(jobs)
        for k,v in ray.get(ready[0]).items():
            if k in kmers:
                kmers[k] += v
            else:
                kmers[k] = v
    if len(kmers):
        print(f"Significant k-mers: {len(kmers)}")
        with open(out_file, 'w') as writer:
            print('k-mer', f'{basename}_Count', sep='\t', file=writer)
            for kmer,count in sorted(kmers.items()):
                print(kmer, count, sep='\t', file=writer)
        return basename,out_file
    else:
        print(f"No significant k-mers found")
        return basename,None

## Create Figures
#
def createFigures(tsv_list:dict, type_string:str, out_path:os.PathLike, lowmem=None, class_file=None, pca=False):
    print(f"\nCreating {type_string} Graphs")

    figPlots = dict()
    start_time = timeit.default_timer()

    # Stacked Bar Plots (top kmer counts)
    combined_tsv = os.path.join(out_path, f'combined_{type_string}.tsv')
    if not os.path.exists(combined_tsv):
        mercat2_report.merge_tsv(tsv_list, combined_tsv)
        if DEBUG:
            print(f"\nTime to merge TSV files: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")

    combined_T_tsv = os.path.join(out_path, f'combined_{type_string}_T.tsv')
    if not os.path.exists(combined_T_tsv):
        mercat2_report.merge_tsv_T(tsv_list, combined_T_tsv)
        if DEBUG:
            print(f"\nTime to merge TSV transposed files: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")


    figPlots[f"Combined {type_string} kmer Summary"] = mercat2_figures.kmer_summary(combined_tsv)
    if DEBUG:
        print(f"Time to compute Combined {type_string} kmer Summary: {round(timeit.default_timer() - start_time,2)} seconds")
        print(f"Virtual Memory {mem_use()}GB")

    # PCA
    if pca and len(tsv_list) > 3:
        print("\nRunning PCA")
        start_time = timeit.default_timer()

        out_pca = os.path.join(out_path, f'pca_{type_string}')
        os.makedirs(out_pca, exist_ok=True)
        pca3d,pca2d = mercat2_figures.plot_PCA(combined_T_tsv, out_pca, lowmem, class_file, DEBUG)
        if pca3d:
            figPlots[f"{type_string} PCA 3D"] = pca3d
        if pca2d:
            figPlots[f"{type_string} PCA 2D"] = pca2d

    return figPlots


## Main
#
def mercat_main():
    '''Main method'''

    # Process Arguments
    args, parser = parseargs()
    # Set Flags
    m_kmer = args.k
    m_num_cores = args.n
    m_inputfile = args.i
    m_inputfolder = args.f
    m_min_count = args.c
    m_flag_prodigal = args.prod
    m_flag_fgs = args.fgs
    m_chunk_size = args.s
    m_outputfolder = Path(args.o)
    m_lowmem = None if args.lowmem is None else bool(args.lowmem)
    m_skipclean = args.skipclean
    m_toupper = args.toupper
    m_class_file = args.category_file
    m_pca = args.pca
    
    if m_outputfolder.exists():
        if args.replace:
            shutil.rmtree(m_outputfolder)
        else:
            parser.error(f"Output folder exists, please specify another folder or use the flag '-replace' to override the files. '{m_outputfolder}'")
    
    m_outputfolder.mkdir(0o777, True, True)

    print(f"\nStarting MerCat2 v{__version__} with k-mer {m_kmer} and {m_num_cores} threads\n")

    ray.init(num_cpus=m_num_cores, log_to_driver=DEBUG)
    if DEBUG:
        print(f"\nVirtual Memory {mem_use()}GB")

    # Load input files
    gc_content = dict()
    samples = dict(nucleotide=dict(),
                protein=dict(),
                prodigal=dict(),
                fgs=dict())

    cleanpath = os.path.join(m_outputfolder, 'clean')
    print("Loading files")
    start_time = timeit.default_timer()
    jobsQC = list()
    jobsFastq = list()
    jobsContig = list()

    @ray.remote(num_cpus=1)
    def fastq_qc(file, cleanpath, basename):
        return (basename, mercat2_fasta.qc(file, cleanpath, basename))
    @ray.remote(num_cpus=1)
    def clean_contig(file, cleanpath, basename, toupper, skipclean):
        if skipclean:
            stat = None
        else:
            file,stat = mercat2_fasta.removeN(file, cleanpath, toupper)
        return (basename, file, stat)
    @ray.remote(num_cpus=1)
    def fastq_to_fasta(file, cleanpath, basename, skiptrim:bool):
        jobsQC = [fastq_qc.remote(file, cleanpath, basename)]
        if not skiptrim:
            file = mercat2_fasta.trim(file, cleanpath, basename)
            jobsQC += [fastq_qc.remote(file, cleanpath, basename)]
        return (basename, mercat2_fasta.fq2fa(file, cleanpath, basename), jobsQC)
    # Load files from folder
    if m_inputfolder:
        m_inputfolder = os.path.abspath(os.path.expanduser(m_inputfolder))
        for fname in os.listdir(m_inputfolder):
            file = Path(m_inputfolder, fname)
            if file.is_file:
                suffixes = file.suffixes
                for i in reversed(range(len(suffixes))):
                    if ''.join(suffixes[i:]) in FILE_EXT_FASTQ + FILE_EXT_NUCLEOTIDE + FILE_EXT_PROTEIN:
                        m_inputfile.append(file)
                        break
    # Read input files
    if m_inputfile:
        for filename in m_inputfile:
            basepath = Path(filename).expanduser().absolute()
            suffixes = basepath.suffixes
            f_ext = ''
            for i in reversed(range(len(suffixes))):
                if ''.join(suffixes[i:]) in FILE_EXT_FASTQ + FILE_EXT_NUCLEOTIDE + FILE_EXT_PROTEIN:
                    f_ext = ''.join(suffixes[i:])
            basename = basepath.name.removesuffix(f_ext)
            print(basename, f_ext)
            if f_ext in FILE_EXT_FASTQ:
                jobsFastq += [fastq_to_fasta.remote(filename, cleanpath, basename, m_skipclean)]
            elif f_ext in FILE_EXT_NUCLEOTIDE:
                jobsContig += [clean_contig.remote(filename, cleanpath, basename, m_toupper, m_skipclean)]
                command = [ 'countAssembly.py', '-f', filename, '-i', '100' ]
                statfile = os.path.join(m_outputfolder, 'stats', f'{basename}.txt')
                os.makedirs(os.path.join(m_outputfolder, 'stats'), exist_ok=True)
                with open(statfile, 'w') as writer:
                    subprocess.run(command, stdout=writer, stderr=subprocess.DEVNULL)
            elif f_ext in FILE_EXT_PROTEIN:
                samples['protein'][basename] = [filename]

    # Wait for jobs
    while jobsFastq:
        ready,jobsFastq = ray.wait(jobsFastq)
        basename,file,jobsqc = ray.get(ready[0])
        samples['nucleotide'][basename] = [file]
        jobsQC += jobsqc
        #gc_content[basename] = stat['GC Content']
    while jobsContig:
        ready,jobsContig = ray.wait(jobsContig)
        basename,file,stat = ray.get(ready[0])
        samples['nucleotide'][basename] = [file]
        if stat:
            gc_content[basename] = stat['GC Content']

    print(f"Time to load {len(samples['nucleotide'])+len(samples['protein'])} files: {round(timeit.default_timer() - start_time,2)} seconds")
    if DEBUG:
        print(f"Virtual Memory {mem_use()}GB")
    start_time = timeit.default_timer()

    # Begin processing files
    figPlots = dict()
    tsv_stats = dict()
    jobsDiversity = list()
    report_dir = os.path.join(m_outputfolder, "report")
    os.makedirs(report_dir, exist_ok=True)

    ## Process Nucleotides ##
    if samples['nucleotide']:
        # Chunk Nucleotides
        if m_chunk_size > 0:
            print("Checking for large nucleotide files")
            start_time = timeit.default_timer()
            dir_chunks = os.path.join(m_outputfolder, "chunks_nucleotide")
            jobs = []
            for basename,file in samples['nucleotide'].items():
                chunk_path = os.path.join(dir_chunks, basename)
                jobs += [chunk_files.remote(basename, *file, m_chunk_size, chunk_path)]
            while jobs:
                ready,jobs = ray.wait(jobs)
                name,chunks = ray.get(ready[0])
                chunks = [Path(c) for c in chunks]
                samples['nucleotide'][name] += chunks
            if DEBUG:
                print(f"Time to check for large files: {round(timeit.default_timer() - start_time,2)} seconds")

        print("Processing Nucleotides")
        print(f"Running Mercat2 using {m_num_cores} cores")
        out_tsv = os.path.join(m_outputfolder, 'tsv_nucleotide')
        os.makedirs(out_tsv, exist_ok=True)
        jobs = []
        start_time = timeit.default_timer()
        for basename,files in samples['nucleotide'].items():
            out_counts = os.path.join(out_tsv, f"{basename}_counts.tsv")
            files = files[1:] if len(files)>1 else [files[0]]
            jobs += [run_mercat2.remote(basename, files, out_counts, m_kmer, m_min_count, m_num_cores)]
        tsv_list = dict()
        while(jobs):
            ready,jobs = ray.wait(jobs)
            basename,kmers = ray.get(ready[0])
            if kmers:
                tsv_list[basename] = kmers
        print(f"Time to count {m_kmer}-mers: {round(timeit.default_timer() - start_time,2)} seconds")
        if DEBUG:
            print(f"Virtual Memory {mem_use()}GB")
        # Stacked Bar Plots (top kmer counts)
        start_time = timeit.default_timer()
        if len(tsv_list):
            figPlots.update(createFigures(tsv_list, "Nucleotide", m_outputfolder, m_lowmem, m_class_file, m_pca))
            tsvfile = Path(m_outputfolder, f'combined_{sample_type}_T.tsv')
            sample_type = "Nucleotide"
            mercat2_diversity.compute_beta_diversity(sample_type, tsvfile, Path(m_outputfolder, "report", "diversity"))
            if DEBUG:
                print(f"Time to calculate Beta Diversity: {round(timeit.default_timer() - start_time,2)} seconds")
        for basename,filename in tsv_list.items():
            os.makedirs(os.path.join(report_dir, 'diversity'), exist_ok=True)
            outfile = os.path.join(report_dir, 'diversity', f'nucleotide-{basename}.tsv')
            jobsDiversity += [diversity.remote(basename, filename, outfile, 'Nucleotide')]
        if len(gc_content) > 0:
            figPlots['Sample GC Summary'] = mercat2_figures.GC_plot_sample(gc_content)


    ## Process Proteins ##
    # PROD ORF Call
    if m_flag_prodigal and samples['nucleotide']:
        @ray.remote(num_cpus=1)
        def orf_call_prod(basename, file, prodpath):
            return mercat2_fasta.orf_call(basename, file, prodpath)

        print(f"\nRunning Prodigal on {len(samples['nucleotide'])} files")
        start_time = timeit.default_timer()
        prodpath = os.path.join(m_outputfolder, 'prodigal')
        os.makedirs(prodpath, exist_ok=True)
        jobsProd = []
        for basename,files in samples['nucleotide'].items():
            jobsProd += [orf_call_prod.remote(basename, files[0], prodpath)]
        while jobsProd:
            ready,jobsProd = ray.wait(jobsProd)
            ret = ray.get(ready[0])
            if ret:
                samples['prod'][ret[0]] = [ret[1]]
        if DEBUG:
            print(f"Time to run Prodigal: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")

    # FGS ORF Call
    if m_flag_fgs and samples['nucleotide']:
        @ray.remote(num_cpus=1)
        def orf_call_fgs(basename, file, outpath):
            return mercat2_fasta.orf_call_fgs(basename, file, outpath)
        print(f"\nRunning FragGeneScanRS on {len(samples['nucleotide'])} files")
        start_time = timeit.default_timer()
        outfgs = os.path.join(m_outputfolder, 'fgs')
        jobsFGS = []
        for basename,files in samples['nucleotide'].items():
            jobsFGS += [orf_call_fgs.remote(basename, files[0], outfgs)]
        while jobsFGS:
            ready,jobsFGS = ray.wait(jobsFGS)
            ret = ray.get(ready[0])
            if ret:
                samples['fgs'][ret[0]] = [ret[1]]
        if DEBUG:
            print(f"Time to run FGS: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")


    ## Process Proteins ##
    for sample_type in samples.keys():
        if (not samples[sample_type]) or (sample_type == "nucleotide"):
            continue
        print(f"\nProcessing Proteins ({sample_type})")
        print(f"Running Mercat2 using {m_num_cores} cores")
        # Chunk large files
        if m_chunk_size > 0:
            print("Checking for large protein files")
            start_time = timeit.default_timer()
            dir_chunks = os.path.join(m_outputfolder, f"chunks_{sample_type}")
            jobs = []
            for basename,file in samples[sample_type].items():
                chunk_path = os.path.join(dir_chunks, basename)
                jobs += [chunk_files.remote(basename, file[0], m_chunk_size, chunk_path)]
            while jobs:
                ready,jobs = ray.wait(jobs)
                name,chunks = ray.get(ready[0])
                chunks = [Path(c) for c in chunks]
                samples[sample_type][name] += chunks
        if DEBUG:
            print(f"Time to check for large protein files: {round(timeit.default_timer() - start_time,2)} seconds")

        jobs = []
        tsv_list = dict()
        out_tsv = os.path.join(m_outputfolder, f"tsv_{sample_type}")
        os.makedirs(out_tsv, exist_ok=True)
        start_time = timeit.default_timer()
        for basename,files in samples[sample_type].items():
            out_counts = os.path.join(out_tsv, f"{basename}_counts.tsv")
            files = files[1:] if len(files)>1 else [files[0]]
            jobs += [run_mercat2.remote(basename, files, out_counts, m_kmer, m_min_count, m_num_cores)]
        while(jobs):
            ready,jobs = ray.wait(jobs)
            if ready[0]:
                basename,kmers = ray.get(ready[0])
                if kmers:
                    tsv_list[basename] = kmers
        print(f"Time to count {m_kmer}-mers: {round(timeit.default_timer() - start_time,2)} seconds")
        if DEBUG:
            print(f"Virtual Memory {mem_use()}GB")
        start_time = timeit.default_timer()
        if len(tsv_list):
            figPlots.update(createFigures(tsv_list, sample_type, m_outputfolder, m_lowmem, m_class_file, m_pca))
            tsvfile = Path(m_outputfolder, f'combined_{sample_type}_T.tsv')
            mercat2_diversity.compute_beta_diversity(sample_type, tsvfile, Path(m_outputfolder, "report", "beta_diversity"))
            if DEBUG:
                print(f"Time to calculate Beta Diversity: {round(timeit.default_timer() - start_time,2)} seconds")
        for basename,filename in tsv_list.items():
            os.makedirs(os.path.join(report_dir, 'diversity'), exist_ok=True)
            outfile = os.path.join(report_dir, 'diversity', f'{sample_type}-{basename}.tsv')
            jobsDiversity += [diversity.remote(basename, filename, outfile, sample_type)]


    # Plot Data
    mercat2_report.write_html(os.path.join(report_dir, "report.html"), figPlots, tsv_stats)
    for sample_type in ['protein', 'fgs', 'prod']:
        if len(samples[sample_type]):
            tsv_out = os.path.join(report_dir, f'metrics-{sample_type}.tsv')
            htm_out = os.path.join(report_dir, f'metrics-{sample_type}.html')
            figPlots = mercat2_figures.plot_sample_metrics(samples[sample_type], tsv_out)
            mercat2_report.write_html(htm_out, figPlots, dict())

    # Wait for any remaining QC Jobs
    if jobsQC:
        print("Waiting for any remaining QC jobs")
    while jobsQC:
        ready,jobsQC = ray.wait(jobsQC)

    print("Gathering Diversity Metrics")
    mergedDiversity = dict()
    while jobsDiversity:
        ready,jobsDiversity = ray.wait(jobsDiversity)
        key,outfile = ray.get(ready[0])
        if key not in mergedDiversity:
            mergedDiversity[key] = []
        mergedDiversity[key].append(outfile)

    for key,filelist in mergedDiversity.items():
        if len(filelist) < 2:
            continue
        tomerge = dict()
        for filename in filelist:
            name = re.search(r'\w+-(\w+).tsv', Path(filename).name)
            if name:
                name = name.group(1)
                tomerge[name] = filename
        outfile = Path(report_dir, f'diversity-{key}.tsv')
        mercat2_report.merge_tsv(tomerge, outfile)

    print("\nFinished MerCat2 Pipeline")

    return


# If main app
if __name__ == "__main__":
    mercat_main()
