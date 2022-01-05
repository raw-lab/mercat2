# -*- coding: utf-8 -*-
"""rheaQC.py: Module for checking quality of .fastq files
Uses FastQC [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/]
$ fastqc file.fastq
$ fastqc file_R1.fastq fastqc file_R2.fastq
"""

import os
import subprocess
import re
import textwrap


## checkQuality
#
def checkQuality(rawRead, config, subdir):
    if type(rawRead) is str:
        return checkSingleRead(rawRead, config, subdir)
    else:
        return checkPairedRead(rawRead, config, subdir)


## checkSingleQuality
#
def checkSingleRead(singleRead, config, subdir):
    path = f"{config['DIR_OUT']}/{subdir}"
    os.makedirs(path, exist_ok=True)
    
    command = f"{config['EXE_FASTQC']} -o {path} {singleRead}"
    try:
        with open(f"{path}/stdout.txt", 'w') as fout, open(f"{path}/stderr.txt", 'w') as ferr:
            subprocess.run(command, shell=True, check=True, stdout=fout, stderr=ferr)
        return os.path.join(path, os.path.splitext(os.path.basename(singleRead))[0]+'_fastqc.html')
    except Exception as e:
        print(e)

    return None


## checkPairedQuality
#
def checkPairedRead(pairedRead, config, subdir):
    path = f"{config['DIR_OUT']}/{subdir}"
    os.makedirs(path, exist_ok=True)
    
    command = f"{config['EXE_FASTQC']} -o {path} {pairedRead[0]} {pairedRead[1]}"
    try:
        with open(f"{path}/stdout.txt", 'w') as fout, open(f"{path}/stderr.txt", 'w') as ferr:
            subprocess.run(command, shell=True, check=True, stdout=fout, stderr=ferr)
    except Exception as e:
        print(e)

    return path

## Split Sequence by N
#
def split_sequenceN(name, sequence):
    N_lengths = []
    regex = re.compile(r"(N+)")
    for match in regex.finditer(sequence):
        N_lengths.append(len(match.group(1)))
    sequences = regex.sub('\n', sequence).split('\n')
    name = name.split()
    basename = name[0]
    info = ' '.join(name[1:])
    seqs = []
    for i, seq in enumerate(sequences, 1):
        header = f">{basename}_{i} {info}"
        seqs.append(header)
        seqs += textwrap.wrap(seq, 80)
    
    return seqs, N_lengths


## Remove N's
#
def removeN(fasta, outpath):
    os.makedirs(outpath, exist_ok=True)    

    outFasta, ext = os.path.splitext(fasta)
    outFasta = os.path.basename(outFasta) + "_clean"+ ext
    outFasta = os.path.join(outpath, outFasta)

    proc = subprocess.run(['grep', '-cE', '^[^>].*N', fasta], stdout=subprocess.PIPE, text=True)
    res = int(proc.stdout.strip())
    if res == 0:
        return fasta, None

    with open(fasta) as reader, open(outFasta, 'w') as writer:
        NStats = dict()
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
                    NStats[name] = stats
                    print('\n'.join(sequences), file=writer)
                else:
                    print('>', name, sep='', file=writer)
                    print('\n'.join(textwrap.wrap(sequence, 80)), file=writer)
                continue #already got next line, next item in loop
            line = reader.readline()

    return os.path.abspath(outFasta), NStats


## Process Fastq
#
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
