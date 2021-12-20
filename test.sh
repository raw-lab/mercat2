#!/bin/bash

#conda_dep = "fastqc fastp prodigal"
#pypi_dep = "setuptools numpy pandas humanize plotly psutil joblib dask scikit-learn scikit-bio"

# source this file to test mercat2
ABSPATH="$(pwd -P)"
export PATH="$ABSPATH/bin:$PATH"
export PYTHONPATH="$ABSPATH:$PYTHONPATH"
