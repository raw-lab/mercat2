#!/bin/bash

#SBATCH --job-name=22-11-09
#SBATCH --partition=Draco
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=0
#SBATCH --time=1-0
#SBATCH -o slurm-%x-%j.out
#SBATCH --mail-type=END,FAIL,REQUEUE

echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""


SECONDS=0

module load anaconda3
eval "$(conda shell.bash hook)"
conda activate mercat2-dev

pip install ~/raw-lab/mercat2

export PYTHONUNBUFFERED=1

echo "##### RUNNING MerCat2 k=4 c=0 #####"
command time mercat2.py -k 4 -f ../data/5_genomes -o $SLURM_JOB_NAME-k4-min0 -fgs -c 0
echo "##### RUNNING MerCat2 k=4 c=3 #####"
command time mercat2.py -k 4 -f ../data/5_genomes -o $SLURM_JOB_NAME-k4-min3 -fgs -c 3
echo "##### RUNNING MerCat2 k=4 c=5 #####"
command time mercat2.py -k 4 -f ../data/5_genomes -o $SLURM_JOB_NAME-k4-min5 -fgs -c 5
echo "##### RUNNING MerCat2 k=4 c=10 ##### "
command time mercat2.py -k 4 -f ../data/5_genomes -o $SLURM_JOB_NAME-k4-min10 -fgs -c 10


echo "##### RUNNING MerCat2 k=31 c=0 #####"
command time mercat2.py -k 31 -f ../data/5_genomes -o $SLURM_JOB_NAME-k31-min0 -fgs -c 0
echo "##### RUNNING MerCat2 k=31 c=3 #####"
command time mercat2.py -k 31 -f ../data/5_genomes -o $SLURM_JOB_NAME-k31-min3 -fgs -c 3
echo "##### RUNNING MerCat2 k=31 c=5 #####"
command time mercat2.py -k 31 -f ../data/5_genomes -o $SLURM_JOB_NAME-k31-min5 -fgs -c 5
echo "##### RUNNING MerCat2 k=31 c=10 #####"
command time mercat2.py -k 31 -f ../data/5_genomes -o $SLURM_JOB_NAME-k31-min10 -fgs -c 10


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "Run Time   : $SECONDS Seconds"
echo "======================================================"
echo ""
