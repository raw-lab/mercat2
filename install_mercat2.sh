#!/usr/bin/env bash
set -e

# initialize conda environment in bash script
eval "$(conda shell.bash hook)"

# create the mercat environment in conda
ENV_NAME=mercat2
conda create -n $ENV_NAME -c conda-forge -c bioconda grpcio python fastqc fastp prodigal metaomestats ray-core ray-dashboard ray-default configargparse pandas numpy humanize plotly psutil dominate scikit-learn scikit-bio scipy==1.8.1 python-kaleido -y
conda activate $ENV_NAME
pip install .
