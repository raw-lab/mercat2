#!/bin/bash

kmer=$1
threads=$2
dataset=$3
runnum=$4
outname="${dataset}-mercat2-${kmer}-${threads}-${runnum}"
outplace="out/outputs/${outname}.out"
tempplace="tmp"

mercat2.py -f data/${dataset}/files -k ${kmer} -n ${threads} -c 10 -o $outplace
