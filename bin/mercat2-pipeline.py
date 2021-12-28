#!/usr/bin/env python

"""mercat2.py: Python code for Parallel k-mer counting."""

__author__      = "Jose L. Figueroa III, Richard A. White III"
__copyright__   = "Copyright 2022"

import sys
import os
import glob
import psutil
import shutil
import timeit
import humanize
import subprocess
import pandas as pd
from collections import OrderedDict
from joblib import Parallel, delayed

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import dask.dataframe as dd

# Mercat libraries
from mercat2 import (mercat2_Chunker, mercat2_report, mercat2_metrics)


# GLOBAL VARIABLES
protein_file_ext = ['.fa','.fna','.ffn','.fasta']


def fastq_processing(fq_path, path, f_name,file):
    trim_path = path+'/'+f_name+"_trim.fastq"
    cmd = "fastqc "+fq_path
    cmd1 = "fastp -i "+fq_path+" -o " + trim_path
    cmd2 = "fastqc "+trim_path
    trim_fna = path+'/'+f_name+"_trim.fna"
    cmd3 = "sed -n '1~4s/^@/>/p;2~4p' "+trim_path+"  > "+trim_fna
    subprocess.call(cmd,shell=True)
    subprocess.call(cmd1,shell=True)
    subprocess.call(cmd2,shell=True)
    subprocess.call(cmd3,shell=True)
    return trim_fna


def check_command(cmd):
    cmd1 = cmd
    with open(os.devnull, 'w') as FNULL:
        try:
            subprocess.check_call(cmd1, stdout=FNULL, stderr=FNULL, shell=True)
        except subprocess.CalledProcessError as e:
            # print e.output -- null since we suppressed output in check_call
            print(("Mercat Error: %s not found, please setup %s using: conda install %s" %(cmd,cmd,cmd)))
            sys.exit(1)


def parseargs(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    num_cores = psutil.cpu_count(logical=False)
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-i', type=str, required=False, help='path-to-input-file') #default=nucleotide
    parser.add_argument('-f', type=str, required=False, help='path-to-folder-containing-input-files')
    parser.add_argument('-k', type=int, required = True, help='kmer length')
    parser.add_argument('-n', type=int, default=num_cores, help='no of cores [default = all]')  # no of cores to use
    parser.add_argument('-c', type=int, default=10, help='minimum kmer count [default = 10]')  # minimum kmer count to report
    parser.add_argument('-pro', action='store_true', help='protein input file')
    parser.add_argument('-p', action='store_true', help='run prodigal on fasta file')
    parser.add_argument('-s', type=int, nargs='?', const=100, required=False, help='Split into x MB files. Default = 100MB')

    # Process arguments
    args = parser.parse_args()
    print(args.i)
    if args.i:
        path_q, file_q = os.path.split(args.i)
        f_name, f_ext = os.path.splitext(file_q)
        if f_ext=='.fastq':
            f_name=f_name.split('.')[0]
            print(path_q+","+args.i)
            args.i=fastq_processing(args.i,path_q,f_name,file_q)
    if args.i and args.f:
        parser.error("Can only specify either an input file (-i) or path to folder containing input files (-f) at a time")
    if args.i:
        if os.path.exists(args.i) < 1:
            path_error = "file " + args.i + " does not exist.\n"
            parser.error(path_error)
    elif args.f:
        if os.path.exists(args.f) < 1:
            path_error = "folder " + args.f + " does not exist.\n"
            parser.error(path_error)
    else:
        parser.error("Please provide either an input file (-i) or an input folder (-f)")
    if args.p:
        if args.pro: parser.error("Can only provide one of -p or -pro option at a time")
        check_command('prodigal')
    return [args,parser]


def calculateKmerCount(cseq,kmer): ###(seq,cseq, prune_kmer,kmer):
    kmerlist = dict()
    for i in range((len(cseq)-kmer)+1):
        k = cseq[i:i+kmer]
        if k not in kmerlist:
            kmerlist[k] = 0
        kmerlist[k] += 1
    return kmerlist


def check_args(ipfile,args,def_option,m_parser):
    given_ext = (os.path.splitext(ipfile)[1]).strip()
    if def_option:
        if given_ext not in protein_file_ext:
            m_parser.error("Input file provided should be one of the following formats: " + str(protein_file_ext))

    if args.p:
        # if not args.q and given_ext not in protein_file_ext:
        if given_ext not in protein_file_ext:
            m_parser.error("Input file provided should be one of the following formats: " + str(protein_file_ext))

    
    if args.pro:
        if given_ext != ".faa":
            m_parser.error("Input file provided should be in .faa format")


def mercat_main():
    __args__, m_parser = parseargs()

    kmer = __args__.k
    num_cores = __args__.n
    m_inputfile = __args__.i
    m_inputfolder = __args__.f
    prune_kmer = __args__.c
    # mflag_fastq = __args__.q
    mflag_prodigal = __args__.p
    # mflag_trimmomatic = __args__.t
    mflag_protein = __args__.pro
    mfile_size_split = __args__.s
    kmerstring = str(kmer) + "-mers"
    
    if not mfile_size_split:
        mfile_size_split = 100

    np_string = "nucleotide"
    if mflag_protein or mflag_prodigal: np_string = "protein"
    # def_option =  not __args__.p and not __args__.q and not __args__.pro
    def_option =  not __args__.p and not __args__.pro

    all_ipfiles = []
    if m_inputfolder:
        m_inputfolder = os.path.abspath(m_inputfolder)
        os.chdir(m_inputfolder)
        #Assume all have same ext
        for fname in os.listdir(m_inputfolder):
            mip = os.path.join(m_inputfolder, fname)
            if not os.path.isdir(mip):
                # skip directories
                path_q, file_q = os.path.split(mip)
                f_name, f_ext = os.path.splitext(file_q)
                if f_ext == '.fastq':
                    f_name = f_name.split('.')[0]
                    # print(path_q+","+args.i)
                    mip = fastq_processing(mip, path_q, f_name, file_q)
                all_ipfiles.append(mip)
    else:
        #m_inputfolder = os.getcwd()
        m_inputfolder = os.path.dirname(os.path.abspath(m_inputfile))
        all_ipfiles.append(os.path.abspath(m_inputfile))
    
    top10_all_samples = dict()
    # print (len(all_ipfiles))
    for m_inputfile in all_ipfiles:
        # if arg
        os.chdir(m_inputfolder)
        check_args(m_inputfile,__args__,def_option,m_parser)

        m_inputfile = os.path.abspath(m_inputfile)

        sample_name = os.path.splitext(os.path.basename(m_inputfile))[0]
        basename_ipfile = os.path.splitext(os.path.basename(m_inputfile))[0] + "_" + np_string

        inputfile_size = os.stat(m_inputfile).st_size
        dir_runs = "mercat_results/" + basename_ipfile + "_run"

        if os.path.exists(dir_runs):
            shutil.rmtree(dir_runs)
        os.makedirs(dir_runs)

        all_chunks_ipfile = []
        is_chunked = False
        if inputfile_size >= (mfile_size_split*1024*1024): #100MB
            print("Large input file provided: Splitting it into smaller files...\n")
            mercat2_Chunker.Chunker(m_inputfile, dir_runs, str(mfile_size_split)+"M", ">")
            os.chdir(dir_runs)
            all_chunks_ipfile = glob.glob("*")
            is_chunked = True
        else:
            os.chdir(dir_runs)
            all_chunks_ipfile.append(m_inputfile)

        splitSummaryFiles = []

        for inputfile in all_chunks_ipfile:

            bif = os.path.splitext(os.path.basename(inputfile))[0] + "_" + np_string

            "Run prodigal if specified"
            '''prodigal -i test_amino-acid.fa -o output.gff -a output.orf_pro.faa  -f gff -p meta -d output.orf_nuc'''
            if mflag_prodigal:
                mflag_protein = True
                gen_protein_file = bif+"_pro.faa"
                prod_cmd = "prodigal -i %s -o %s -a %s -f gff -p meta -d %s" % (
                inputfile, bif + ".gff", gen_protein_file, bif + "_nuc.ffn")
                print(prod_cmd)
                with open(os.devnull, 'w') as FNULL:
                    subprocess.call(prod_cmd, stdout=FNULL, stderr=FNULL, shell=True)
                inputfile = gen_protein_file

            print(("Running mercat using " + str(num_cores) + " cores"))
            print(("input file: " + inputfile))

            start_time = timeit.default_timer()
            sequences = OrderedDict()
            is_fastq = False
            with open(inputfile,'r') as f:
                for line in f:
                    if line.startswith(">"): break
                    elif line.startswith("@"):
                        is_fastq = True
                        break
            with open(inputfile,'r') as f:
                if not is_fastq:
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
                    #assert sname and seq
                    sequences[sname] = seq
                else: #process fastq file
                    seq = ""
                    sname = ""
                    for line in f:
                        line = line.strip()
                        if line.startswith("@"):
                            sname = line[1:].split()[0]
                        elif line.startswith("+"):
                            if seq:
                                sequences[sname] = seq
                                seq = ""
                        else:
                            if sname not in sequences: seq = line

            #print sequences.keys()[0] + "="+ sequences.values()[0]

            print(("Number of sequences in " + inputfile + " = "+ str(humanize.intword(len(sequences)))))
            
            results = Parallel(n_jobs=num_cores)(
                delayed(calculateKmerCount)(sequences[seq], kmer) for seq in sequences)

            kmerlist = dict()
            #kmerlist_all_seq = dict()

            for d in results:
                for k,v in list(d.items()):
                    if k in kmerlist:
                        kmerlist[k] += v
                    else: kmerlist[k] = v

            print(("Time to compute " + kmerstring +  ": " + str(round(timeit.default_timer() - start_time,2)) + " secs"))

            significant_kmers = []
            for k in kmerlist:
                if kmerlist[k] >= prune_kmer:
                    significant_kmers.append(k)

            print(("Total number of " + kmerstring +  " found: " + str(humanize.intword(len(kmerlist)))))
            print((kmerstring +  " with count >= " + str(prune_kmer) + ": " + str(humanize.intword(len(significant_kmers)))))

            if mflag_protein:
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
                    # df.set_value(k, 'Count', kmerlist[k])
                    df.at[k,'Count'] = kmerlist[k]
                    len_cseq = float(len(c_kmer))
                    df.at[k,'GC_Percent'] = round(((c_kmer.count("G")+c_kmer.count("C")) / len_cseq) * 100.0)
                    df.at[k,'AT_Percent'] = round(((c_kmer.count("A")+c_kmer.count("T")) / len_cseq) * 100.0)
            #df.to_csv(bif + "_summary.csv", index_label=kmerstring, index=True)
            df.to_csv(bif + "_summary.csv", index_label="k-mers", index=True)

            splitSummaryFiles.append(bif + "_summary.csv")
            print(("Total time: " + str(round(timeit.default_timer() - start_time,2)) + " secs"))

        num_chunks = len(all_chunks_ipfile)
        df = dd.read_csv(splitSummaryFiles)
        dfgb = df.groupby("k-mers").sum()
        df10 = dfgb.nlargest(10,'Count').compute()
        dfsum = dfgb.sum(0).compute()

        dfgb.to_csv("./" + basename_ipfile + "_finalSummary*.csv", index_label=kmerstring, name_function=lambda name: str(name))

        if mflag_protein:
            df10[['PI', 'MW', 'Hydro']] = df10[['PI', 'MW', 'Hydro']] / num_chunks
        else:
            df10[['GC_Percent', 'AT_Percent']] = df10[['GC_Percent', 'AT_Percent']] / num_chunks

        top10_all_samples[sample_name] = [df10,dfsum.Count]

        all_counts = dfgb.Count.values.compute().astype(int)
        mercat2_metrics.compute_alpha_beta_diversity(all_counts,basename_ipfile)

        if is_chunked:
            for tempfile in all_chunks_ipfile:
                os.remove(tempfile)
            for sf in splitSummaryFiles:
                os.remove(sf)

    plots_dir = os.path.join(m_inputfolder, "mercat_results", "plots")
    if os.path.exists(plots_dir):
        shutil.rmtree(plots_dir)
    os.makedirs(plots_dir)
    os.chdir(plots_dir)

    figPlots = dict()
    tsv_stats = dict()
    # Stacked Bar Plots (top kmer counts)
    if len(all_ipfiles) == 1:
        sbname = os.path.basename(all_ipfiles[0])
    else:
        sbname = os.path.basename(m_inputfolder)
    figPlots[sbname+'_k-mers'] = mercat2_report.stackedbar_plots(sbname, top10_all_samples, 'Count', kmerstring)
    
    # PCA
    subdir = m_inputfolder+"/mercat_results/"
    if len(all_ipfiles) >= 4:
       figPlots['PCA'] = mercat2_report.PCA_plot(subdir)

    # Scatter Plots
    for basename_ipfile in top10_all_samples:
        df10,_ = top10_all_samples[basename_ipfile]
        if mflag_protein:
            figPlots[basename_ipfile+'_PI'] = mercat2_report.scatter_plots(basename_ipfile, 'PI', df10, kmerstring)
            figPlots[basename_ipfile+'_MW'] = mercat2_report.scatter_plots(basename_ipfile, 'MW', df10, kmerstring)
            figPlots[basename_ipfile+'_Hydro'] = mercat2_report.scatter_plots(basename_ipfile, 'Hydro', df10, kmerstring)
        else:
            figPlots[basename_ipfile+'_GC'] = mercat2_report.scatter_plots(basename_ipfile, 'GC_Percent', df10, kmerstring)
            figPlots[basename_ipfile+'_AT'] = mercat2_report.scatter_plots(basename_ipfile, 'AT_Percent', df10, kmerstring)

    # Save HTML Report
    mercat2_report.write_html(os.path.join(plots_dir, "report.html"), figPlots, tsv_stats)
    return 0


if __name__ == "__main__":
    mercat_main()
