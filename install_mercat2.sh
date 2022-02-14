#!/usr/bin/env bash

ABSPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

function install {
    install_path=$1
    mkdir -p $install_path
    cp bin/*.py $install_path
    cp bin/*.sh $install_path
    # TODO: Copy packages as well.
    echo "Program files copied to '$install_path'"
    echo "Add this to your PATH or .bashrc for easier use:"
    echo "export PATH=\"$install_path:\$PATH\""
    return
}

function install_pip {
    rm -r dist/
    echo "Building Mercat2 distribution..."
    python -m build > /dev/null #2>&1
    rm -r *.egg-info/
    # install latest build version
    latest=$(ls dist/*.whl | sort -V | tail -n 1)
    python -m pip uninstall mercat2 -y
    echo
    echo "Installing $latest"
    echo
    python -m pip install $latest
    return
}

function install_conda {
    # initialize conda environment in bash script
    eval "$(conda shell.bash hook)"

    # create the mercat environment in conda
    conda env remove --name mercat -y
    conda create -n mercat -c conda-forge -c bioconda fastqc fastp prodigal python=3.9 -y

    # install additional pip requirements
    conda activate mercat
    pip install build setuptools
    # TODO: Change this to proper install once uploaded to pypi and bioconda.
    install_pip
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
    -i|--install)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        ARG_INSTALL=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -p|--pip)
      ARG_PIP=true
      shift
      ;;
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

    -i, --install PATH  Copy scripts to PATH and downloads the database.
                        Use this option for a manual install.
                        Assumes dependencies are already installed.
                        (requires 'unzip', 'wget', and 'git')
    -p, --pip           Instal Mercat2 using pip from local folder.
                        Dependencies must be installed manually.
                        (requires 'pip')
    -c, --conda         Creates a conda environment named 'Mercat2' with all dependencies and installs Mercat2 in it
                        (requires Anaconda3 to be installed)
    -h, --help          Display this message and exit
" && exit 0

[ $ARG_INSTALL ] && install $ARG_INSTALL && exit 0
[ $ARG_PIP ] && install_pip && exit 0
[ $ARG_CONDA ] && install_conda && exit 0
[ $ARG_ENV ] && develop_env
