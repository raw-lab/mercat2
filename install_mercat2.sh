#!/usr/bin/env bash
set -e

# initialize conda environment in bash script
eval "$(conda shell.bash hook)"

# create the mercat environment in conda
ENV_NAME=mercat2
mamba create -n $ENV_NAME -c bioconda -c conda-forge grpcio'=1.43' python">=3.9" fastqc fastp prodigal metaomestats ray-core ray-dashboard ray-default ray-tune configargparse pandas numpy humanize plotly psutil dominate scikit-learn scikit-bio scipy python-kaleido -y
conda activate $ENV_NAME
pip install .
