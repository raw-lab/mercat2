#!/bin/bash

kmer=$1
threads=$2
dataset=$3
runnum=$4
outname="${dataset}-kmc-${kmer}-${threads}-${runnum}"
outplace="out/outputs/${outname}.out"
tempplace="tmp"

./progs/kmc -b -k${kmer} -ci10 -t${threads} -fm @data/${dataset}/inputfiles.lst ${outplace} ${tempplace}
