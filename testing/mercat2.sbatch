#!/bin/bash

#SBATCH --partition=Andromeda
#SBATCH --job-name=mercat-test
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256GB
#SBATCH --time=1-0
#SBATCH -e slurm-%x-%j.err
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

module load anaconda3
eval "$(conda shell.bash hook)"
conda activate mercat

command time mercat2-pipeline.py -k 4 -f 5_genomes

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
