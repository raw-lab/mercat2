#!/bin/bash

kmer=$1
threads=$2
dataset=$3
runnum=$4
outname="${dataset}-jellyfish-${kmer}-${threads}-${runnum}"
outplace="out/outputs/${outname}.out"
tempplace="tmp"

./progs/jellyfish count -m $kmer -t $threads -L 10 -s 1M -o $outplace data/${dataset}/files/*
