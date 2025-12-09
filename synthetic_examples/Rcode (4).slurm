#!/bin/bash
#SBATCH --job-name=KS1
#SBATCH --array=88
#SBATCH --time=00-01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=16
#SBATCH -p hns

ml R/4.2.0
Rscript syn1.R $SLURM_ARRAY_TASK_ID
