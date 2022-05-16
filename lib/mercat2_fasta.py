# -*- coding: utf-8 -*-
"""mercat2_fasta.py: Module with tools for fasta files
"""

import os
import subprocess
import re
import textwrap


## Split Sequence by N
def split_sequenceN(header:str, sequence:str):
    '''Splits a nucleotide sequence string at occurrences of N.
    Helper method for removeN.
    Assumes the first word before the first space is the name of the sequence.

    Parameters:
        header (str): The header of the sequence, without the leading '>'.
        sequence (str): The nucleotide sequence.

    Returns:
        tuple: A tuple with a list of sequences (including modified headers)
               and a corresponding list of contiguous N repeats in each sequence.
    '''

    N_lengths = []
    regex = re.compile(r"(N+)")
    for match in regex.finditer(sequence):
        N_lengths.append(len(match.group(1)))
    sequences = regex.sub('\n', sequence).split('\n')
    header = header.split()
    basename = header[0]
    info = ' '.join(header[1:])
    seqs = []
    for i, seq in enumerate(sequences, 1):
        header = f">{basename}_{i} {info}"
        seqs.append(header)
        seqs += textwrap.wrap(seq, 80)
    
    return (seqs, N_lengths)


## Remove N's
def removeN(fasta:str, outpath:str):
    '''Splits sequences in a scaffold fasta file at N repeats.

    Parameters:
        fasta (str): A path to a fasta file.
        outpath (str): A folder path of where to save the modified file.

    Returns:
        tuple: A tuple with the path to the cleaned file and a dictionary containing basic stats.
    '''

    os.makedirs(outpath, exist_ok=True)    

    outFasta, ext = os.path.splitext(fasta)
    outFasta = os.path.basename(outFasta) + "_clean"+ ext
    outFasta = os.path.join(outpath, outFasta)

    NStats = dict()
    gc_count = 0
    total_length = 0
    with open(fasta) as reader, open(outFasta, 'w') as writer:
        line = reader.readline()
        while line:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:]
                sequence = ""
                line = reader.readline()
                while line:
                    line = line.strip()
                    if line.startswith('>'):
                        break
                    sequence += line
                    line = reader.readline()
                if 'N' in sequence:
                    sequences, stats = split_sequenceN(name, sequence)
                    print('\n'.join(sequences), file=writer)
                    for seq in sequences:
                        gc_count += seq.count('G') + seq.count('C')
                        total_length += len(seq)
                else:
                    print('>', name, sep='', file=writer)
                    print('\n'.join(textwrap.wrap(sequence, 80)), file=writer)
                    gc_count += sequence.count('G') + sequence.count('C')
                    total_length += len(sequence)
                continue #already got next line, next item in loop
            line = reader.readline()
    NStats['GC Content'] = 100.0 * gc_count / total_length

    return (os.path.abspath(outFasta), NStats)


## Process Fastq
def fastq_processing(fq_file:str, outpath:str, f_name:str):
    '''Processes fastq files.
    Uses fastqc to get the quality of the reads.
    Uses fastp to trim the reads.
    Uses Linux sed command to convert the fastq file to fasta format.

    Parameters:
        fq_file (str): The path to a fastq file.
        outpath (str): The path to save the fastqc report, trimmed fastq file, and fasta file.
        f_name (str): The name of the sample.

    Returns:
        str: The absolute path to the final trimmed fasta file.
    '''

    os.makedirs(outpath, exist_ok=True)
    trim_fq = os.path.join(outpath, f_name+"_trim.fastq")
    trim_fna = os.path.join(outpath, f_name+"_trim.fna")
    
    cmd_qc = f"fastqc {fq_file} -o {outpath}"
    cmd_fastp = f"fastp -i {fq_file} -o {trim_fq}"
    cmd_qc2 = f"fastqc {trim_fq} -o {outpath}"
    cmd_format = f"sed -n '1~4s/^@/>/p;2~4p' {trim_fq} > {trim_fna}"
    subprocess.call(cmd_qc, shell=True)
    subprocess.call(cmd_fastp, shell=True)
    subprocess.call(cmd_qc2, shell=True)
    subprocess.call(cmd_format, shell=True)
    return os.path.abspath(trim_fna)


## ORF Call Prodigal
def orf_call(basename:str, file:str, outpath:str):
    '''Finds the ORFs in the nucleotides and converts to an amino acid file.
    Uses prodigal for ORF calling.
    Produces a protein ffa file.
    Produces a gff file.

    Parameters:
        basename (str): The base name of the file for output files.
        file (str): The path to a nucleotide fasta file.
        outpath (str): The path to save the output files.

    Returns:
        tuple: A tuple with the name, and path to the protein faa file.
    '''

    outpath = os.path.abspath(outpath)
    out_pro = os.path.join(outpath, basename+"_pro.faa")
    out_gff = os.path.join(outpath, basename+".gff")
    out_nuc = os.path.join(outpath, basename+"_nuc.fna")
    prod_cmd = f"prodigal -i {file} -o {out_gff} -a {out_pro} -f gff -p meta"
    os.makedirs(outpath, exist_ok=True)
    with open(os.devnull, 'w') as FNULL:
        subprocess.call(prod_cmd, stdout=FNULL, stderr=FNULL, shell=True)
    return (basename, out_pro)
