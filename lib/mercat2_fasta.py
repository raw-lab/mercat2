# -*- coding: utf-8 -*-
"""mercat2_fasta.py: Module with tools for fasta files
"""

import os
import sys
import pkg_resources as pkg
import shutil
import subprocess
import tarfile
import re
import textwrap

PATH_FGS = pkg.resource_filename("mercat2_lib", "FGS")


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
    with open(fasta, 'r') as reader, open(outFasta, 'w') as writer:
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


## Check Command
#
def check_command(cmd:str):
    '''Checks if the command is available in PATH and executable'''

    if shutil.which(cmd) is None:
        print(f"Mercat Error: {cmd} not found, please setup {cmd} using: 'conda install {cmd}'")
        return False
    else:
        return True


## Process Fastq
def qc(fq_file:str, outpath:str, f_name:str):
    '''Processes fastq files.
    Uses fastqc to get the quality of the reads.

    Parameters:
        fq_file (str): The path to a fastq file.
        outpath (str): The path to save the fastqc report.
        f_name (str): The name of the sample.

    Returns:
        None
    '''

    os.makedirs(outpath, exist_ok=True)
    if check_command('fastqc'):
        subprocess.run(['fastqc', fq_file, '-o', outpath], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return

def trim(fq_file:str, outpath:str, f_name:str):
    '''Processes fastq files.
    Uses fastp to trim the reads.

    Parameters:
        fq_file (str): The path to a fastq file.
        outpath (str): The path to save the trimmed fastq file.
        f_name (str): The name of the sample.

    Returns:
        str: The path to the trimmed fasta file.
    '''

    os.makedirs(outpath, exist_ok=True)
    trim_fq = os.path.join(outpath, f_name+"_trim.fastq")
    if check_command('fastp'):
        subprocess.run(['fastp', '-i', fq_file, '-o', trim_fq], stdout=open(f'{outpath}/{f_name}-trim.stdout', 'w'), stderr=open(f'{outpath}/{f_name}-trim.stderr', 'w'))
    return trim_fq

def fq2fa(fq_file:str, outpath:str, f_name:str):
    '''Processes fastq files.
    Uses Linux sed command to convert the fastq file to fasta format.

    Parameters:
        fq_file (str): The path to a fastq file.
        outpath (str): The path to save the  fasta file.
        f_name (str): The name of the sample.

    Returns:
        str: The path to the converted fasta file.
    '''

    os.makedirs(outpath, exist_ok=True)
    fna_file = os.path.join(outpath, f_name+".fna")

    # convert fastq to fasta
    subprocess.run(['sed', '-n', '1~4s/^@/>/p;2~4p', fq_file], stdout=open(fna_file, 'w'))
    return os.path.abspath(fna_file)


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

    if not check_command('prodigal'):
        exit()

    outpath = os.path.abspath(outpath)
    out_pro = os.path.join(outpath, basename+"_pro.faa")
    #out_gff = os.path.join(outpath, basename+".gff")
    #out_nuc = os.path.join(outpath, basename+"_nuc.fna")
    #prod_cmd = f"prodigal -i {file} -o {out_gff} -a {out_pro} -f gff -p meta"
    prod_cmd = ['prodigal',
                '-i', file,
                '-a', out_pro,
                #'-o', out_gff,
                #'-f', 'gff',
                '-p', 'meta']
    os.makedirs(outpath, exist_ok=True)
    with open(f'{outpath}/{basename}.stdout', 'w') as stdout, open(f'{outpath}/{basename}.stderr', 'w') as stderr:
        subprocess.run(prod_cmd, stdout=stdout, stderr=stderr)
    return out_pro


## ORF Call FragGeneScanRS
def orf_call_fgs(basename:str, file:str, outpath:str):
    '''Finds the ORFs in the nucleotides and converts to an amino acid file.
    Uses FragGeneScanRS for ORF calling.
    Produces a protein ffa file.
    Produces a gff file.

    Parameters:
        basename (str): The base name of the file for output files.
        file (str): The path to a nucleotide fasta file.
        outpath (str): The path to save the output files.

    Returns:
        tuple: A tuple with the name, and path to the protein faa file.
    '''

    exe_fgs = os.path.join(PATH_FGS, 'FragGeneScanRs')

    if not os.path.exists(exe_fgs):
        system = sys.platform

        if system == "Windows":
            print("Windows is not supported with FragGeneScanRs")
            return None

        with tarfile.open(os.path.join(PATH_FGS, f'FragGeneScanRS-{system}.tar.gz'), 'r') as fgs:
            fgs.extractall(PATH_FGS)

    outpath = os.path.abspath(outpath)
    os.makedirs(outpath, exist_ok=True)
    faa_out = os.path.join(outpath, f'{basename}.faa')

    command = ['FragGeneScanRs',
                '--complete',
                '-s', file,
                '-t', 'complete',
                '-a', faa_out,
                ]

    subprocess.run(command)

    return (basename, f"{faa_out}")
