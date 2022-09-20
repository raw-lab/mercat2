#!/usr/bin/env bash
set -e

function install_conda {
    # initialize conda environment in bash script
    eval "$(conda shell.bash hook)"

    # create the mercat environment in conda
    ENV_NAME=mercat2-dev
    mamba create -n $ENV_NAME -c conda-forge -c bioconda fastqc fastp prodigal metaomestats ray-core ray-dashboard ray-default configargparse pandas dask dask-ml numpy humanize plotly psutil dominate scikit-learn python-kaleido -y
    conda activate $ENV_NAME
    pip install .
    return
}

### Begin Main Script ###

install_conda
