#!/usr/bin/env bash

set -e

function install_conda {
    # initialize conda environment in bash script
    eval "$(conda shell.bash hook)"

    # create the mercat environment in conda
    conda create -n mercat2 -c conda-forge -c bioconda fastqc fastp prodigal ray-core ray-dashboard configargparse pandas numpy humanize plotly psutil joblib dominate scikit-learn -y
    conda activate mercat2
    pip install .
    return
}

### Begin Main Script ###

install_conda
