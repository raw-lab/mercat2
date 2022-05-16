#!/usr/bin/env bash

set -e

ABSPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

function install_conda {
    # initialize conda environment in bash script
    eval "$(conda shell.bash hook)"

    # create the mercat environment in conda
    conda create -n mercat2 -c conda-forge -c bioconda fastqc fastp ray-core ray-dashboard prodigal -y

    # install additional pip requirements
    conda activate mercat2
    #pip install build setuptools
    # TODO: Change this to proper install once uploaded to pypi and bioconda.
    pip install .
    return
}

function develop_env {
  # source ./install_mercat2.sh -e
  ABSPATH=$(pwd)
  export PATH="$ABSPATH/bin:$PATH"
  export PYTHONPATH="$ABSPATH:$PYTHONPATH"
  return
}

### Begin Main Script ###

# Parse Arguments
while (( "$#" )); do
  case "$1" in
    -c|--conda)
      ARG_CONDA=true
      shift
      ;;
    -e|--env)
      ARG_ENV=true
      shift
      ;;
    -h|--help)
      ARG_HELP=true
      shift
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done

[ ! $ARG_INSTALL $ARG_PIP $ARG_CONDA $ARG_ENV $ARG_HELP ] && echo "
No options given.
" && ARG_HELP=true

[ $ARG_HELP ] && echo "
usage: [--path PATH] [--download] [--dependencies] [--help]

    -c, --conda         Creates a conda environment named 'Mercat2' with all dependencies and installs Mercat2 in it
                        (requires Anaconda3 to be installed)
    -h, --help          Display this message and exit
" && exit 0

[ $ARG_PIP ] && install_pip && exit 0
[ $ARG_CONDA ] && install_conda && exit 0
[ $ARG_ENV ] && develop_env
