#!/usr/bin/env bash
set -e

# initialize conda environment in bash script
eval "$(conda shell.bash hook)"

# create the mercat environment in conda
ENV_NAME=MerCat2

mamba create -n $ENV_NAME -y -c conda-forge -c bioconda python">=3.9" \
	fastqc fastp metaomestats \
	hydrampp configargparse humanize plotly psutil dominate scikit-learn scikit-bio scipy python-kaleido matplotlib pyrodigal

conda activate $ENV_NAME
pip install .
