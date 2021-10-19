#!/bin/sh
#SBATCH --job-name=decline_random
#SBATCH --partition=batch
#SBATCH --time=10:00:00
module load R/3.5.1
module load slim

Rscript  checkFile3.R

slim random_decline.slim
wait
