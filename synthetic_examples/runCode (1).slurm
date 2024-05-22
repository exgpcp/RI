#!/bin/bash
#SBATCH --job-name=MAPf1
#SBATCH --array=1-300
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G
#SBATCH -p hns

ml R/4.2.0
Rscript weighted MAP for hyperparameters.R $SLURM_ARRAY_TASK_ID
