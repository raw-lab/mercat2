#!/usr/bin/env bash
set -e

# initialize conda environment in bash script
eval "$(conda shell.bash hook)"

# create the mercat environment in conda
ENV_NAME=mercat2-dev
mamba create -n $ENV_NAME -c conda-forge -c bioconda grpcio python'>=3.10' fastqc fastp prodigal metaomestats ray-core ray-dashboard ray-default configargparse pandas dask dask-ml numpy humanize plotly psutil dominate scikit-learn scikit-bio scipy==1.8.1 python-kaleido -y
conda activate $ENV_NAME
pip install .
