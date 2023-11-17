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
progs = [ "mercat2", "kmc", "jellyfish" ]
data_sources = [ "5genome-fna", "archaeal-viruses-82-fna", "cpr-78-fna", "gtdb-archaea-100-fna", "gtdb-bacteria-100-fna", "phages-100-fna", "viruses-100-fna" ]
x = np.arange(3)

data_headers = [
    "Time (Seconds)",
    "RAM (M)",
    "Disk (M)"
]

data_deviations = [
    "time-dev",
    "ram-dev",
    "disk-dev"
]

kmers = [ 4, 31 ]

fancy_title_map = {
    "5genome-fna": "Five Genome",
    "archaeal-viruses-82-fna": "Archaeal Viruses",
    "cpr-78-fna": "CPR",
    "gtdb-archaea-100-fna": "GTDB Archaea",
    "gtdb-bacteria-100-fna": "GTDB Bacteria",
    "phages-100-fna": "Phages",
    "viruses-100-fna": "Viruses",
}


for source in data_sources:
    df = outer_df \
        .filter((pl.col("Dataset") == source)) \
        .group_by("Program", "Kmer Length", maintain_order=True) \
        .agg([
            pl.col("Time (Seconds)"),
            pl.col("time-dev"),
            pl.col("RAM (M)"),
            pl.col("ram-dev"),
            pl.col("Disk (M)"),
            pl.col("disk-dev")])

    plt.rcParams["figure.dpi"] = 120
    fig, ax = plt.subplots(2, 3)
    fig.suptitle(f"{fancy_title_map[source]}", y=0.95, color="#333333", weight="bold")
    fig.set_figwidth(10)
    fig.set_figheight(8)
    fig.tight_layout(pad=3.0)
    fig.patch.set_linewidth(2)
    fig.patch.set_edgecolor("#777777")
    fig.subplots_adjust(hspace=0.5)

    for k in range(2):
        kmer = kmers[k]
        for j in range(3):
            newdf = df.filter(pl.col("Kmer Length") == kmer)
            my_data = newdf[data_headers[j]].to_list()
            my_devs = newdf[data_deviations[j]].to_list()

            if j == 0:
                ax[k][j].set_title(f"k = {kmer}", loc="left", weight="bold", x=-0.1)
            ax[k][j].spines["top"].set_visible(False)
            ax[k][j].spines["right"].set_visible(False)
            ax[k][j].spines["left"].set_visible(False)
            ax[k][j].spines["bottom"].set_color("#DDDDDD")
            ax[k][j].tick_params(bottom=False, left=False)
            ax[k][j].set_axisbelow(True)
            ax[k][j].yaxis.grid(True, color="#EEEEEE")

            ax[k][j].set_xticks(x, threads)
            ax[k][j].set_xlabel("Threads", labelpad=3, color="#333333")
            ax[k][j].set_ylabel(data_headers[j], labelpad=3, color="#333333")
            ax[k][j].bar(x-0.2, my_data[0], 0.2, yerr=my_devs[0], color="#00BDAE", capsize=3)
            ax[k][j].bar(x    , my_data[1], 0.2, yerr=my_devs[1], color="#404041", capsize=3)
            ax[k][j].bar(x+0.2, my_data[2], 0.2, yerr=my_devs[2], color="#357CD2", capsize=3)

    fig.legend(progs, loc="outside upper center", ncol=3, bbox_to_anchor=(0.5, 0.48))
    plt.savefig(f"images/{fancy_title_map[source]}.png", edgecolor=fig.get_edgecolor())
    plt.close()
