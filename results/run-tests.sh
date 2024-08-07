#!/bin/bash

set -ex

PRE=$(date +%Y-%m-%d)
STDOUT=$PRE/stdout.txt
STDERR=$PRE/stderr.txt

rm -rf $PRE
mkdir -p $PRE

PYTHONUNBUFFERED=1

for K in 4 15 31 59
do
	for C in 0 1 2 10
	do
		for S in 1 10
		do
			command time mercat2.py -k $K -c $C -f ../data/5-genomes-fna/ -pca -prod -fgs -s $S -o $PRE/fna-5genomes-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -f ../data/5-genomes-fna_gz/ -pca -prod -fgs -s $S -o $PRE/fna-5genomes_gz-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -f ../data/5-genomes-faa/ -pca -prod -fgs -s $S -o $PRE/faa-5genomes-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -f ../data/5-genomes-faa_gz/ -pca -prod -fgs -s $S -o $PRE/faa-5genomes_gz-k$K-c$C-$S

			command time mercat2.py -k $K -c $C -i ../data/5-genomes-fna/DJ.fna -prod -fgs -s $S -o $PRE/fna-DJ-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -i ../data/5-genomes-fna_gz/DJ.fna.gz -prod -fgs -s $S -o $PRE/fna-DJ_gz-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -i ../data/5-genomes-faa/DJ_pro.faa -prod -fgs -s $S -o $PRE/faa-DJ_gz-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -i ../data/5-genomes-faa_gz/DJ_pro.faa.gz -prod -fgs -s $S -o $PRE/DJ_gz-k$K-c$C-$S

			command time mercat2.py -k $K -c $C -i ../data/5-genomes-fna/RW* -pca -prod -fgs -s $S -o $PRE/fna-5genomes-i-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -i ../data/5-genomes-fna_gz/RW* -pca -prod -fgs -s $S -o $PRE/fna-5genomes_gz-i-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -i ../data/5-genomes-faa/RW* -pca -prod -fgs -s $S -o $PRE/faa-5genomes-i-k$K-c$C-$S
			command time mercat2.py -k $K -c $C -i ../data/5-genomes-faa_gz/RW* -pca -prod -fgs -s $S -o $PRE/faa-5genomes_gz-i-k$K-c$C-$S
		done

		command time mercat2.py -k $K -c $C -i ../data/Test_R1.fastq -prod -fgs -o $PRE/test-qc-k$K-c$C
		command time mercat2.py -k $K -c $C -i ../data/Test_R1.fastq.gz -prod -fgs -o $PRE/test-qc_gz-k$K-c$C
	done
done
