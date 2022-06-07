#!/bin/bash

#SBATCH --job-name=22-06-03_prodigal_5-genomes
#SBATCH --partition=Draco
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=6
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
conda activate mercat2-ray

command time mercat2.py -prod -k 4 -f ../data/5_genomes -o "results_$SLURM_JOB_NAME"


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "Run Time   : $SECONDS Seconds"
echo "======================================================"
echo ""
