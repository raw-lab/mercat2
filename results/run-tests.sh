#!/bin/bash

set -x

PRE=$(date +%Y-%m-%d)
STDOUT=$PRE/stdout.txt
STDERR=$PRE/stderr.txt

rm -r $PRE
mkdir -p $PRE

PYTHONUNBUFFERED=1

for chunk in 1 10
do
	command time mercat2.py -k 5 -f ../data/5-genomes-fna/ -pca -prod -fgs -s $chunk -o $PRE/fna-5genomes-$chunk
	command time mercat2.py -k 5 -f ../data/5-genomes-fna_gz/ -pca -prod -fgs -s $chunk -o $PRE/fna-5genomes_gz-$chunk
	command time mercat2.py -k 5 -f ../data/5-genomes-faa/ -pca -prod -fgs -s $chunk -o $PRE/faa-5genomes-$chunk
	command time mercat2.py -k 5 -f ../data/5-genomes-faa_gz/ -pca -prod -fgs -s $chunk -o $PRE/faa-5genomes_gz-$chunk

	command time mercat2.py -k 5 -i ../data/5-genomes-fna/DJ.fna -prod -fgs -s $chunk -o $PRE/fna-DJ-$chunk
	command time mercat2.py -k 5 -i ../data/5-genomes-fna_gz/DJ.fna.gz -prod -fgs -s $chunk -o $PRE/fna-DJ_gz-$chunk
	command time mercat2.py -k 5 -i ../data/5-genomes-faa/DJ_pro.faa -prod -fgs -s $chunk -o $PRE/faa-DJ_gz-$chunk
	command time mercat2.py -k 5 -i ../data/5-genomes-faa_gz/DJ_pro.faa.gz -prod -fgs -s $chunk -o $PRE/DJ_gz-$chunk
done

command time mercat2.py -k 5 -i ../data/Test_R1.fastq -prod -fgs -o $PRE/test-qc
command time mercat2.py -k 5 -i ../data/Test_R1.fastq.gz -prod -fgs -o $PRE/test-qc_gz
