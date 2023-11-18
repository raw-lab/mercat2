#!/usr/bin/env python3
import subprocess
import re
from typing import Tuple
import statistics
import os
script_dir = os.path.dirname(__file__)

kmers = [ 4, 31 ]
datasets = [ "5genome-fna", "archaeal-viruses-82-fna", "cpr-78-fna", "gtdb-archaea-100-fna", "gtdb-bacteria-100-fna", "phages-100-fna", "viruses-100-fna" ]
programs = [ "mercat2", "kmc", "jellyfish" ]
threads = [ 1, 4, 8 ]
runs = [ 1, 2, 3, 4, 5 ]

# Speed > brevity :)
timereg = re.compile(r' seconds$')
timereg2 = re.compile(r's$')
timereg3 = re.compile(r'^.*0:')
ramreg = re.compile(r'^.*: ')
diskreg = re.compile(r' *total.*$')
diskreg2 = re.compile(r'MB$')

outlines = []

def mercat2_fun(outpath: str, termpath: str, timepath: str) -> Tuple[float, float, float]:
    time_output = ram_output = disk_output = -1.0
    with open(termpath) as termfile:
        for line in termfile:
            if "Time to count" in line:
                t = ramreg.sub('', line.strip())
                t = timereg.sub('', t)
                time_output = float(t)
                break
    with open(timepath, 'r') as timefile:
        for line in timefile:
            if "Maximum resident set size" in line:
                r = ramreg.sub('', line.strip())
                ram_output = float(r) / 1000
                break
    disk_output = float(subprocess.check_output(['du','-s', str(outpath)]).split()[0].decode('utf-8')) / 1000

    return (time_output, ram_output, disk_output)

def kmc_fun(outpath: str, termpath: str, timepath: str) -> Tuple[float, float, float]:
    time_output = ram_output = disk_output = -1.0
    with open(termpath) as termfile:
        for line in termfile:
            if "Total    :" in line:
                t = ramreg.sub('', line.strip())
                t = timereg2.sub('', t)
                time_output = float(t)
                continue
            if "Tmp size" in line:
                ds = ramreg.sub('', line.strip())
                ds = diskreg2.sub('', ds)
                disk_output = float(ds)
                break
    with open(timepath) as timefile:
        for line in timefile:
            if "Maximum resident set size" in line:
                r = ramreg.sub('', line.strip())
                ram_output = float(r) / 1000
                break

    return (time_output, ram_output, disk_output)

def jellyfish_fun(outpath: str, termpath: str, timepath: str) -> Tuple[float, float, float]:
    time_output = ram_output = disk_output = -1.0
    with open(timepath, 'r') as timefile:
        for line in timefile:
            if "Elapsed (wall clock)" in line:
                tmp = ramreg.sub('', line.strip())
                tmp = tmp.split(':')
                mins = int(tmp[0])
                secs = float(tmp[1])
                t = mins * 60 + float(secs)
                time_output = t
                continue
            if "Maximum resident set size" in line:
                r = ramreg.sub('', line.strip())
                ram_output = float(r) / 1000
                break

    disk_output = float(subprocess.check_output(['du','-s', str(outpath)]).split()[0].decode('utf-8')) / 1000

    return (time_output, ram_output, disk_output)

func_map = {
        "mercat2": mercat2_fun,
        "kmc": kmc_fun,
        "jellyfish": jellyfish_fun,
}

for dataset in datasets:
    for program in programs:
        for kmer in kmers:
            for thread in threads:
                times = []
                ramsizes = []
                disksizes = []
                for runnum in runs:
                    outname = f'{dataset}-{program}-{kmer}-{thread}-{runnum}'
                    outpath = os.path.join(script_dir, "out/outputs/", f"{outname}.out")
                    timepath = os.path.join(script_dir, "out/timeinfo/", f"{outname}.timeinfo")
                    termpath = os.path.join(script_dir, "out/termout/", f"{outname}.termout")

                    time_output, ram_output, disk_output = func_map[program](outpath, termpath, timepath)
                    times.append(time_output)
                    ramsizes.append(ram_output)
                    disksizes.append(disk_output)

                time_avg = statistics.fmean(times)
                time_dev = statistics.pstdev(times)
                ram_avg = statistics.fmean(ramsizes)
                ram_dev = statistics.pstdev(ramsizes)
                disk_avg = statistics.fmean(disksizes)
                disk_dev = statistics.pstdev(disksizes)
                print(dataset, program, kmer, thread, times, ramsizes, disksizes)
                outlines.append(f"{dataset},{kmer},{program},{thread},{time_avg},{time_dev},{ram_avg},{ram_dev},{disk_avg},{disk_dev}\n")

column_names = "Dataset,Kmer Length,Program,Threads,Time (Seconds),time-dev,RAM (M),ram-dev,Disk (M),disk-dev\n"
with open("results.csv", 'w') as writefile:
    writefile.write(column_names)
    writefile.writelines(outlines)

