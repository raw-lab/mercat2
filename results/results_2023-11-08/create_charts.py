#!/usr/bin/env python3

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg") # so saving images works

import argparse
parser = argparse.ArgumentParser()
required = parser.add_argument_group('Required arguments')
required.add_argument('-i', '--input', type=str, help='Path to input csv', required=True)

args = parser.parse_args()

outer_df = pl.read_csv(args.input)
threads = [ '1', '4', '8' ]
kmers = [ 4, 31 ]
progs = [ "mercat2", "kmc", "jellyfish" ]
data_sources = [ "5genome-fna", "archaeal-viruses-82-fna", "cpr-78-fna", "gtdb-archaea-100-fna", "gtdb-bacteria-100-fna", "phages-100-fna", "viruses-100-fna" ]
x = np.arange(3)

plt.rcParams["figure.dpi"] = 120

for source in data_sources:
    for kmer in kmers:
        df = outer_df \
            .filter(
                    (pl.col("Kmer Length") == kmer) &
                    (pl.col("Dataset") == source)) \
            .group_by("Program", maintain_order=True) \
            .agg([
                pl.col("Time (Seconds)"),
                pl.col("time-dev"),
                pl.col("RAM (M)"),
                pl.col("ram-dev"),
                pl.col("Disk (M)"),
                pl.col("disk-dev")])

        data = [
            "Time (Seconds)",
            "RAM (M)",
            "Disk (M)"
        ]

        devs = [
            "time-dev",
            "ram-dev",
            "disk-dev"
        ]

        fancy_title_map = {
            "5genome-fna": "5 Genome",
            "archaeal-viruses-82-fna": "Archaeal Viruses 82 Genome",
            "cpr-78-fna": "CPR 78 Genome",
            "gtdb-archaea-100-fna": "GTDB Archaea 100 Genome",
            "gtdb-bacteria-100-fna": "GTDB Bacteria 100 Genome",
            "phages-100-fna": "Phages 100 Genome",
            "viruses-100-fna": "Viruses 100 Genome",
        }

        fig, ax = plt.subplots(1, 3)
        fig.suptitle(f"{fancy_title_map[source]} K={kmer}", y=0.95, color="#333333", weight="bold")
        fig.set_figwidth(10)
        fig.tight_layout(pad=3.0)
        fig.patch.set_linewidth(2)
        fig.patch.set_edgecolor("#777777")

        for i in range(3):
            my_data = df[data[i]].to_list()
            my_devs = df[devs[i]].to_list()

            ax[i].spines["top"].set_visible(False)
            ax[i].spines["right"].set_visible(False)
            ax[i].spines["left"].set_visible(False)
            ax[i].spines["bottom"].set_color("#DDDDDD")
            ax[i].tick_params(bottom=False, left=False)
            ax[i].set_axisbelow(True)
            ax[i].yaxis.grid(True, color="#EEEEEE")

            ax[i].set_xticks(x, threads)
            ax[i].set_xlabel("Threads", labelpad=3, color="#333333")
            ax[i].set_ylabel(data[i], labelpad=3, color="#333333")
            ax[i].bar(x-0.2, my_data[0], 0.2, yerr=my_devs[0], color="#00BDAE", capsize=3)
            ax[i].bar(x, my_data[1], 0.2, yerr=my_devs[1], color="#404041", capsize=3)
            ax[i].bar(x+0.2, my_data[2], 0.2, yerr=my_devs[2], color="#357CD2", capsize=3)

        fig.legend(progs, loc="outside upper center", ncol=3, bbox_to_anchor=(0.5, 0.9))
        plt.savefig(f"images/{fancy_title_map[source]}-{kmer}.png", edgecolor=fig.get_edgecolor())
        plt.close()
