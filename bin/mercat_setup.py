import io
import re
import subprocess
import sys
import os
import getpass

import shutil
import tempfile
from os.path import isfile, join

############################Set paths for file interactions###################################

path = "/home/%s/Desktop/mercat2"% getpass.getuser()
# path2 = "/home/%s/Desktop/mercat2/osf_Files"% getpass.getuser()
path3 = "/home/%s/Desktop/gittemp"% getpass.getuser()
path_to_wrapper = "/home/%s/Desktop/gittemp/bin/"% getpass.getuser()
access_rights = 0o755

############################Creates the mercat2 folder########################################

def mercat2_dir():
    try:
        os.mkdir(path, access_rights)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print("Successfully created the directory %s" % path),
        # osf_Files_dir()

if __name__ == "__mercat2_dir__":
    mercat2_dir()


mercat2_dir()

################################Install dependencies#########################

def install_dependencies():

    env_cmd = "conda create -n mercat_env -c conda-forge -c bioconda scikit-bio dask setuptools pandas numpy humanize plotly psutil joblib prodigal scikit-learn"
    subprocess.call(env_cmd, shell=True)
    
    
    git_cmd = "pip install GitPython"
    subprocess.call(git_cmd, shell=True)

    fastqc_cmd = "sudo apt install fastqc"
    subprocess.call(fastqc_cmd, shell=True)

    fastp_cmd = "conda install -c bioconda fastp"
    subprocess.call(fastp_cmd, shell=True)

if __name__ == "__install_dependencies__":
    install_dependencies()

install_dependencies()
#############################get current wrapper from github###################################
import git
def wrapper_download():

    git_URL = 'https://github.com/raw-lab/mercat2.git'
    os.mkdir(path3, access_rights)
    os.chdir(path3)
    git.Repo.clone_from(git_URL, path3, branch='master')
    shutil.move(os.path.join(path_to_wrapper, 'mercat'), path)
    shutil.rmtree(path3)

if __name__ == "__wrapper_download__":
    wrapper_download()

wrapper_download()

