#!/bin/sh
#SBATCH --job-name=nWF_random
#SBATCH --partition=batch
#SBATCH --time=10:00:00
module load R/3.5.1
module load slim

Rscript  checkFile.R

slim random_nWF_nocouples.slim
wait
