#!/bin/bash
#SBATCH --job-name=sny1BM
#SBATCH --array=1-300
#SBATCH --time=0-01:00:00
#SBATCH --mem=6G
#SBATCH -p stat

ml python/3.9.0
python3 BMsyn1_2.py $SLURM_ARRAY_TASK_ID
