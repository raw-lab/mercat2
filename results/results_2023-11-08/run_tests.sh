#!/bin/bash

clear_cache() {
	sudo swapoff -a
	sleep 5
	sudo swapon -a
	sync && sleep 5 && sync
	sudo sysctl -w vm.drop_caches=3 > /dev/null
	sleep 10
}

stop_cron() {
	sudo systemctl stop crond.service
}

start_cron() {
	sudo systemctl start crond.service
}

clean_tmp() {
	rm -rf tmp/*
}

arrThreads=(1 4 8)
arrProgs=(mercat2 jellyfish kmc)
arrK=(4 31)

datasets=(5genome-fna archaeal-viruses-82-fna cpr-78-fna gtdb-archaea-100-fna gtdb-bacteria-100-fna phages-100-fna viruses-100-fna)

eval "$(conda shell.bash hook)"
conda activate mercat2

echo "Stopping cron..."
stop_cron

for dataset in ${datasets[@]} ; do
for kmer in ${arrK[@]} ; do
for prog in ${arrProgs[@]} ; do
for thread in ${arrThreads[@]} ; do
for runnum in {1..5} ; do
	printf "Clearing cache... "
	clean_tmp
	clear_cache
	echo "Done!"
	printf "Running ${prog} with ${dataset} for length ${kmer} with ${thread} threads: trial ${runnum}..."
	command time -v ./helper/${prog}.sh $kmer $thread $dataset $runnum > out/termout/${dataset}-${prog}-${kmer}-${thread}-${runnum}.termout 2> out/timeinfo/${dataset}-${prog}-${kmer}-${thread}-${runnum}.timeinfo
	echo "Done!"
done
done
done
done
done

echo "Starting cron again..."
start_cron
